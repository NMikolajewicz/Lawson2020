

import pandas as pd
import numpy as np
from scipy.stats import ttest_1samp
from statsmodels.stats.multitest import multipletests


def run_multcompare(df,pval_name):
    """
    wrapper for multicomparison correction (adds fdr to a df with a list of pvals)
    """

    pvals = np.array(df[pval_name].values);
    z = pvals.shape
    pvals = pvals.flatten();
    idx = ~np.isnan(pvals)
    idx = np.where(idx)
    ms = multipletests(pvals[idx],alpha=0.05,method='fdr_bh')
    qq = np.ones(len(pvals))
    rej = np.array([False for x in pvals])
    qq[idx] =ms[1]
    rej[idx] = ms[0]
    df['fdr'] = qq
    df['reject'] = rej
    return df

# load data
drug_data = pd.read_csv('TNF_autophinib_growthRates.csv')
# use fraction of decrease in growth rate as metric
drug_data['1 - mu/mu_max'] = 1 - drug_data['mu/mu_max']

synergy_score = [];
cell_lines = [];
exp_ids = [];
dfs = []
synergy_err = []
tstats = []
pvals = []


# for each experiment, compute the synergy score using the bliss independence model

for exp_id,df in drug_data.groupby('exp_id'):
    # for each experiment, build dataframes for  1-mu/mu_max for single agent treatments, the 
    cell_line = df['cell_line'].iloc[0]
    auto_only = df[df['TNF alpha'] == 0]
    auto_only = auto_only[['Autophinib','1 - mu/mu_max']]
    auto_only.set_index(['Autophinib'],inplace=True)
    auto_only.columns = ['1 - mu/mu_max (AUT)']

    tnf_only = df[df['Autophinib'] == 0]
    tnf_only = tnf_only[['TNF alpha','1 - mu/mu_max']]
    tnf_only.set_index(['TNF alpha'],inplace=True)
    tnf_only.columns = ['1 - mu/mu_max (TNF)']

    errs = []
    cpds = []
    for cpd,dfr in tnf_only.groupby('TNF alpha'):
        err = np.abs(dfr['1 - mu/mu_max (TNF)'].iloc[1] - dfr['1 - mu/mu_max (TNF)'].iloc[0])
        errs.append(err)
        cpds.append(cpd)

    tnf_err = pd.DataFrame({'TNF alpha': cpds, 'TNF err': errs})

    errs = []
    cpds = []
    for cpd,dfr in auto_only.groupby('Autophinib'):
        err = np.abs(dfr['1 - mu/mu_max (AUT)'].iloc[1] - dfr['1 - mu/mu_max (AUT)'].iloc[0])
        errs.append(err)
        cpds.append(cpd)

    aut_err = pd.DataFrame({'Autophinib': cpds, 'AUT err': errs})
    dff = df[(df['TNF alpha'] > 0) & (df['Autophinib'] > 0)]
    
    
    # join datafame with combo treat with single agent arms, and compute residual between measured growth rate decrease, and the the sume of single agent 1-mu/mu_max
    df_final = dff.set_index('TNF alpha').join(tnf_only).reset_index().set_index('Autophinib').join(auto_only).reset_index()
    df_final = df_final.set_index('TNF alpha').join(tnf_err.set_index('TNF alpha')).reset_index()
    df_final = df_final.set_index('Autophinib').join(aut_err.set_index('Autophinib')).reset_index()
    df_final['error'] = np.sqrt(df_final['TNF err']**2 + df_final['AUT err']**2)
    # compute the residual between the growth rate decrease using both compounds simutaneously vs. the growth rate decrease with the combination treatment
    df_final['residual'] = df_final['1 - mu/mu_max'] - (df_final['1 - mu/mu_max (TNF)'] + df_final['1 - mu/mu_max (AUT)'])

    # compute synergy scores
    synergy = df_final.residual.mean()
    synergy_std = df_final.residual.std()
    residuals = df_final.residual
    # perform a t-test to see whether residual distribution deviates from population with zero mean
    tstat,pval = ttest_1samp(residuals,0)    
    dfs.append(df_final)
    exp_ids.append(exp_id)
    cell_lines.append(cell_line)
    synergy_score.append(synergy)
    synergy_err.append(synergy_std)
    tstats.append(tstat)
    pvals.append(pval)

res = pd.DataFrame({'exp_id':exp_ids,'Cell Line':cell_lines,'Synergy Score (mean)':synergy_score,'Synergy Score (std)':synergy_err, 'T-statistic': tstats, 'p-value': pvals}).set_index('exp_id')
# run multiple comparison
res = run_multcompare(res,'p-value')
res.to_csv('SynergyScores.csv')
