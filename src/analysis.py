"""Run with local private inputs; see docs/input_format.md."""
from pathlib import Path
import argparse,json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from models import (fit_lorentzian_amplitude, parabolic_peak_refine,
                    estimate_background, find_half_power_points,
                    weighted_linear_fit, interpolate_width_uncertainty,
                    exp_decay_fit_with_offset)
ROOT=Path(__file__).resolve().parents[1]

def run(root=ROOT):
    root=Path(root);folder=root/'private_data'
    names=['settings.json','cw_scans.csv','pulses.csv']
    if any(not (folder/name).is_file() for name in names):
        print('Private measurements are not distributed. Complete the templates in private_data/.')
        return None
    cfg=json.loads((folder/'settings.json').read_text())
    required=['cw_method','L_cw_m','u_L_cw_m','L_pulse_m','u_L_pulse_m','u_delay_us',
              'excluded_modes','frequency_width_fraction','frequency_error_floor_kHz',
              'linewidth_fractional_error','pulse_decay_fractional_error']
    if any(cfg.get(k) is None for k in required):raise ValueError('Complete required settings; see docs/input_format.md.')
    if cfg['cw_method'] not in ['lorentzian','crossings']:raise ValueError('Choose lorentzian or crossings.')
    for k in ['L_cw_m','L_pulse_m','u_delay_us','frequency_width_fraction','frequency_error_floor_kHz']:
        if not np.isfinite(cfg[k]) or cfg[k]<=0:raise ValueError(k+' must be positive and finite.')
    for k in ['u_L_cw_m','u_L_pulse_m','linewidth_fractional_error','pulse_decay_fractional_error']:
        if not np.isfinite(cfg[k]) or cfg[k]<0:raise ValueError(k+' must be nonnegative and finite.')
    if not isinstance(cfg['excluded_modes'],list):raise ValueError('excluded_modes must be a list.')
    scans=pd.read_csv(folder/'cw_scans.csv');pulses=pd.read_csv(folder/'pulses.csv').sort_values('pulse_number')
    for df in [scans,pulses]:
        if df.empty or not np.isfinite(df.to_numpy(dtype=float)).all():raise ValueError('Supply finite numeric measurements.')
    if (scans[['mode','frequency_kHz']]<=0).any().any():raise ValueError('Mode and frequency must be positive.')
    if (pulses.amplitude_mV<=0).any():raise ValueError('Pulse amplitudes must be positive.')
    if len(pulses)<3 or pulses.pulse_number.duplicated().any() or not (np.diff(pulses.delay_us)>0).all():raise ValueError('Supply at least three distinct successive pulses with increasing delays.')
    upto=cfg.get('fit_upto')
    if upto is not None and (not isinstance(upto,int) or not 3<=upto<=len(pulses)):raise ValueError('fit_upto must be null or a valid count of at least three pulses.')
    rows=[]
    for mode,g in scans.groupby('mode',sort=True):
        g=g.sort_values('frequency_kHz');f=g.frequency_kHz.to_numpy();v=g.amplitude_mV.to_numpy()
        if len(f)<5 or len(np.unique(f))!=len(f) or np.ptp(v)<=0:raise ValueError('Each scan needs at least five unique frequencies and varying amplitude.')
        peak=int(np.argmax(v))
        if peak in [0,len(v)-1]:raise ValueError('Scan must bracket its peak.')
        if cfg['cw_method']=='lorentzian':
            r=fit_lorentzian_amplitude(f,v)
            center,width,width_error=r['f0_kHz'],r['delta_f_kHz'],2*r['u_gamma_kHz']
        else:
            center,_=parabolic_peak_refine(f,v)
            _,left,right=find_half_power_points(f,v,estimate_background(v),peak)
            if left is None or right is None:raise ValueError('Half-power crossings not covered.')
            width=right-left;width_error=np.nan
        rows.append(dict(mode=mode,center_kHz=center,width_full_kHz=width,width_fit_error_kHz=width_error,Q=center/width))
    cw=pd.DataFrame(rows)
    extra=folder/'observed_modes_optional.csv'
    if extra.is_file():
        observed=pd.read_csv(extra).sort_values('mode')
        if observed.empty or not np.isfinite(observed.to_numpy(dtype=float)).all():raise ValueError('Remove empty optional file or supply mode measurements.')
        n=observed['mode'].to_numpy();freq=observed.frequency_kHz.to_numpy()
        _,sigma=interpolate_width_uncertainty(cw['mode'].to_numpy(),cw.width_full_kHz.to_numpy(),n)
    else:
        n=cw['mode'].to_numpy();freq=cw.center_kHz.to_numpy()
        sigma=np.maximum(cfg['frequency_width_fraction']*cw.width_full_kHz.to_numpy(),cfg['frequency_error_floor_kHz'])
    mask=~np.isin(n,cfg['excluded_modes'])
    if mask.sum()<3 or len(np.unique(n[mask]))<3:raise ValueError('At least three distinct modes must remain in the speed fit.')
    m,b,dm,db,_,_,chi=weighted_linear_fit(n[mask],freq[mask],sigma[mask])
    if m<=0:raise ValueError('CW slope must be positive.')
    c_cw=2*cfg['L_cw_m']*1000*m
    uc=c_cw*np.hypot(dm/m,cfg['u_L_cw_m']/cfg['L_cw_m'])
    # Preserve original calculation explicitly as legacy: full/half-width needs review.
    alpha=4*np.pi*cw.width_full_kHz.to_numpy()*1000/c_cw
    we=np.fmax(cw.width_fit_error_kHz.to_numpy(),cfg['linewidth_fractional_error']*cw.width_full_kHz.to_numpy())
    cw['alpha_cw_legacy_per_m']=alpha
    cw['u_alpha_cw_legacy_per_m']=alpha*np.hypot(we/cw.width_full_kHz.to_numpy(),uc/c_cw)
    num=pulses.pulse_number.to_numpy();t=pulses.delay_us.to_numpy();amp=pulses.amplitude_mV.to_numpy()
    mt,bt,dmt,_,yt,rt,chit=weighted_linear_fit(num,t,np.full_like(t,cfg['u_delay_us'],dtype=float))
    if mt<=0:raise ValueError('Pulse time slope must be positive.')
    c_p=2*cfg['L_pulse_m']/(mt*1e-6)
    up=c_p*np.hypot(dmt/mt,cfg['u_L_pulse_m']/cfg['L_pulse_m'])
    offset,a0,rate,model,r2=exp_decay_fit_with_offset(t,amp,fit_upto=upto)
    ap=2*rate*1e6/c_p;uap=ap*np.hypot(cfg['pulse_decay_fractional_error'],up/c_p)
    selected=~cw['mode'].isin(cfg['excluded_modes'])
    summary=dict(c_cw_m_per_s=c_cw,u_c_cw=uc,c_pulse_m_per_s=c_p,u_c_pulse=up,
                 cw_reduced_chi2=chi,pulse_reduced_chi2=chit,
                 pulse_decay_rate_per_us=rate,pulse_offset_mV=offset,pulse_amplitude0_mV=a0,
                 pulse_decay_r2_all_points=r2,alpha_pulse_per_m=ap,u_alpha_pulse=uap,
                 alpha_cw_legacy_mean_per_m=float(np.mean(alpha)),
                 alpha_cw_legacy_selected_mean_per_m=float(np.mean(alpha[selected])),
                 velocity_difference_m_per_s=c_p-c_cw,
                 velocity_difference_over_independent_errors=(c_p-c_cw)/np.hypot(uc,up),
                 dispersion_residual_indicator=float(np.max(np.abs(freq-(m*n+b)))/np.mean(sigma)))
    out=root/'private_outputs';out.mkdir(exist_ok=True)
    cw.to_csv(out/'cw_parameters.csv',index=False)
    (out/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False))
    fig,ax=plt.subplots(2,2,figsize=(11,8))
    for mode,g in scans.groupby('mode'):ax[0,0].plot(g.frequency_kHz,g.amplitude_mV,'o-',label=f'mode {mode:g}')
    ax[0,0].set(xlabel='Frequency (kHz)',ylabel='Amplitude (mV)',title='CW scans');ax[0,0].legend()
    ax[0,1].errorbar(n,freq,yerr=sigma,fmt='o');ax[0,1].plot(n,m*n+b)
    ax[0,1].set(xlabel='Mode number',ylabel='Frequency (kHz)',title='CW speed fit')
    ax[1,0].errorbar(num,t,yerr=cfg['u_delay_us'],fmt='o');ax[1,0].plot(num,yt)
    ax[1,0].set(xlabel='Pulse number',ylabel='Delay (microseconds)',title='Pulse speed fit')
    ax[1,1].plot(t,amp,'o',label='Measured');ax[1,1].plot(t,model,label='Offset exponential')
    ax[1,1].set(xlabel='Delay (microseconds)',ylabel='Amplitude (mV)',title='Pulse decay');ax[1,1].legend()
    fig.tight_layout();fig.savefig(out/'analysis.png',dpi=160);plt.close(fig)
    print('Analysis completed. Local results: private_outputs/. CW attenuation uses the original legacy width convention; see methods notes.')
    return summary,cw

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--root',type=Path,default=ROOT)
    try:run(parser.parse_args().root)
    except (ValueError,KeyError,RuntimeError) as error:parser.exit(1,str(error)+'\n')
