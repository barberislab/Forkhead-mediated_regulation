import pandas as pd
import sqlite3
import numpy as np
from scipy.signal import find_peaks

###################################################
# CONNECT TO THE GEMMER DATABASE
###################################################
def create_connection(db_file):
    """ create a database connection to the SQLite database
        specified by db_file
    :param db_file: database file
    :return: Connection object or None
    """
    try:
        conn = sqlite3.connect(db_file)
        return conn
    except Exception as e:
        print((e.message, e.args))
 
    return None

# connect to database
database = "../Data/DB_genes_and_interactions_2022-02-27.db"
conn = create_connection(database)

print('=== Connection established ===')

query = "SELECT standard_name, systematic_name FROM genes"
df_conv = pd.read_sql_query(query, conn).set_index('standard_name')

query = "SELECT standard_name, systematic_name, desc, expression_peak_phase, expression_peak_time, is_enzyme, KEGG_pathway FROM genes"
df_genes = pd.read_sql_query(query, conn).set_index('systematic_name')


###################################################
# INITIATE ALL RELEVANT COLUMNS
###################################################
df_genes['Fkh1 MacIsaac'] = False
df_genes['Fkh1 Venters'] = False
df_genes['Fkh1 Ostrow'] = False
df_genes['Fkh1 Mondeel'] = False
df_genes['Fkh2 MacIsaac'] = False
df_genes['Fkh2 Venters'] = False
df_genes['Fkh2 Ostrow'] = False
df_genes['Fkh2 Mondeel'] = False


###################################################
# GET TARGET STATUS FROM GEMMER FOR: Mondeel et al., Ostrow et al. and MacIsaac et al.
###################################################
query = "SELECT * \
        FROM interactions \
        WHERE type = 'regulation' \
        AND source = 'FKH1'"
df_fkh1 = pd.read_sql_query(query, conn)

query = "SELECT * \
        FROM interactions \
        WHERE type = 'regulation' \
        AND source = 'FKH2'"
df_fkh2 = pd.read_sql_query(query, conn)

MacIsaac = df_conv.loc[df_fkh1.loc[df_fkh1['evidence'].str.contains('16522208'),'target'],'systematic_name']
Ostrow = df_conv.loc[df_fkh1.loc[df_fkh1['evidence'].str.contains('24504085'),'target'],'systematic_name']
Mondeel = df_conv.loc[df_fkh1.loc[df_fkh1['evidence'].str.contains('31299083'),'target'],'systematic_name']

df_genes.loc[MacIsaac,'Fkh1 MacIsaac'] = True
df_genes.loc[Ostrow,'Fkh1 Ostrow'] = True
df_genes.loc[Mondeel,'Fkh1 Mondeel'] = True

MacIsaac = df_conv.loc[df_fkh2.loc[df_fkh2['evidence'].str.contains('16522208'),'target'],'systematic_name']
Ostrow = df_conv.loc[df_fkh2.loc[df_fkh2['evidence'].str.contains('24504085'),'target'],'systematic_name']
Mondeel = df_conv.loc[df_fkh2.loc[df_fkh2['evidence'].str.contains('31299083'),'target'],'systematic_name']

df_genes.loc[MacIsaac,'Fkh2 MacIsaac'] = True
df_genes.loc[Ostrow,'Fkh2 Ostrow'] = True
df_genes.loc[Mondeel,'Fkh2 Mondeel'] = True

###################################################
# VENTERS ET AL.
###################################################
# Because SGD reports the heat shock activated genes (which do not necessarily show up in UTMAX 25C) there are too many targets from GEMMER
# restrict for this analysis to just significant genes from UTmax 25C

df_venters = pd.read_excel('../Data/Venters_25C_UTmax.xls', header=0, skiprows=[0,1,2,4,5,6,7,8,9,10,11,12], usecols=[0,8,9])
df_venters = df_venters.set_index('Factor') # pandas automatically labels the first column with systematic names 'Factor'

venters_fkh1 = df_venters[df_venters['Fkh1'] > 0.8].index.tolist()
venters_fkh2 = df_venters[df_venters['Fkh2'] > 1.13].index.tolist()

venters_fkh1 = [g for g in venters_fkh1 if g in df_genes.index.tolist()]
venters_fkh2 = [g for g in venters_fkh2 if g in df_genes.index.tolist()]

df_genes.loc[venters_fkh1,'Fkh1 Venters'] = True
df_genes.loc[venters_fkh2,'Fkh2 Venters'] = True

###################################################
# ROSSI ET AL.
###################################################
df_genes['Fkh1 Rossi'] = False
df_genes['Fkh2 Rossi'] = False

df_rossi = pd.read_excel('../Data/Rossi2021.xlsx', engine='openpyxl', index_col=0)

Rossi  = [g for g in df_rossi.query('Fkh1 == 1').index.tolist() if g in df_genes.index.tolist()]
df_genes.loc[Rossi,'Fkh1 Rossi'] = True

Rossi  = [g for g in df_rossi.query('Fkh2 == 1').index.tolist() if g in df_genes.index.tolist()]
df_genes.loc[Rossi,'Fkh2 Rossi'] = True


###################################################
# LUPO ET AL.
###################################################
df_lupo = pd.read_csv('../Data/Lupo2021_SumProm_data.csv', header=0, index_col=1)
df_lupo = df_lupo[['cer Fkh1', 'cer Fkh2']] 

df_lupo['Fkh1 Lupo z-score'] = (df_lupo['cer Fkh1'] - df_lupo['cer Fkh1'].mean()) / df_lupo['cer Fkh1'].std()
df_lupo['Fkh2 Lupo z-score'] = (df_lupo['cer Fkh2'] - df_lupo['cer Fkh2'].mean()) / df_lupo['cer Fkh2'].std()
df_lupo = df_lupo.drop(['cer Fkh1', 'cer Fkh2'],axis=1)

df_genes = df_genes.merge(df_lupo, how='left', left_index=True, right_index=True, sort=True)

# define targets as having a z-score > 0
df_genes['Fkh1 Lupo'] = df_genes['Fkh1 Lupo z-score'] > 0
df_genes['Fkh2 Lupo'] = df_genes['Fkh2 Lupo z-score'] > 0


###################################################
# COUNT NUMBER OF TIMES LISTED AS TARGET IN THE 6 STUDIES
###################################################
df_genes.loc[:,'Fkh1_times_target'] = df_genes.loc[:,['Fkh1 MacIsaac','Fkh1 Venters', 'Fkh1 Ostrow','Fkh1 Mondeel','Fkh1 Rossi','Fkh1 Lupo']].sum(axis=1)
df_genes.loc[:,'Fkh2_times_target'] = df_genes.loc[:,['Fkh2 MacIsaac','Fkh2 Venters', 'Fkh2 Ostrow','Fkh2 Mondeel','Fkh2 Rossi','Fkh2 Lupo']].sum(axis=1)


###################################################
# HACKETT ET AL.
###################################################
df_hackett = pd.read_excel('../Data/Hackett2020.xlsx', engine='openpyxl')
df_hackett.columns = ['Gene','Fkh1_OE_0', 'Fkh1_OE_5', 'Fkh1_OE_8', 'Fkh1_OE_15', 'Fkh1_OE_30', 'Fkh1_OE_45', 'Fkh1_OE_60', 'Fkh1_OE_90',
                      'Fkh2_OE_0', 'Fkh2_OE_5', 'Fkh2_OE_8', 'Fkh2_OE_15', 'Fkh2_OE_30', 'Fkh2_OE_45', 'Fkh2_OE_60', 'Fkh2_OE_90']

# deduce if the timecourse is consistently up or down
def effect_fkh_overexpression(row):

    if (row == 0).all():
        return 'Invariant'
    elif row.abs().max() < np.log2(1.1): # max. less than 10% up or down
        return 'Weak response'

    # disregard entries with less than 10% change 
    # i.e. don't consider time points with log2(0.9) < FC < log2(1.1)
    # this catches silly cases where there is one time point going against the trend that is the opposite sign but small
    filtered_row = row[(row.abs() >= np.log2(1.1))]
    
    if (filtered_row >= 0).all(): # upward trend
        if filtered_row.max() >= 1:
            verdict = 'Strongly up'
        else:
            verdict = 'Up'
    elif (filtered_row <= 0).all(): # downward trend
        if filtered_row.min() <= -1:
            verdict = 'Strongly down'
        else:
            verdict = 'Down'
    else: # both up and down by at least 10%
        # try to infer order of peaks
        min_idx = np.argmin(filtered_row)
        max_idx = np.argmax(filtered_row)
        min_logFC = filtered_row.min()
        max_logFC = filtered_row.max()
        if min_idx < max_idx:
            if min_logFC < -1 and max_logFC > 1:
                verdict = 'Strongly down then strongly up'
            elif min_logFC < -1 and max_logFC <= 1:
                verdict = 'Strongly down then up'
            elif min_logFC >= -1 and max_logFC > 1:
                verdict = 'Down then strongly up'
            else:
                verdict = 'Down then up'
        else:
            if min_logFC < -1 and max_logFC > 1:
                verdict = 'Strongly up then strongly down'
            elif min_logFC < -1 and max_logFC <= 1:
                verdict = 'Up then strongly down'
            elif min_logFC >= -1 and max_logFC > 1:
                verdict = 'Strongly up then down'
            else:
                verdict = 'Up then down'
    return verdict
    
df_hackett['Fkh1_OE'] = df_hackett[['Fkh1_OE_0', 'Fkh1_OE_5', 'Fkh1_OE_8', 'Fkh1_OE_15', 'Fkh1_OE_30', 'Fkh1_OE_45', 'Fkh1_OE_60', 'Fkh1_OE_90']].apply(effect_fkh_overexpression,axis=1)
df_hackett['Fkh2_OE'] = df_hackett[['Fkh2_OE_0', 'Fkh2_OE_5', 'Fkh2_OE_8', 'Fkh2_OE_15', 'Fkh2_OE_30', 'Fkh2_OE_45', 'Fkh2_OE_60', 'Fkh2_OE_90']].apply(effect_fkh_overexpression,axis=1)

# some entries in the 'Gene' column are systematic names like YBR139W (ATG42). Try to convert as many as possible to their standard name
replacement_series = df_hackett[df_hackett['Gene'].isin(df_conv['systematic_name'])][['Gene']]
replacement_series['standard_name'] = replacement_series['Gene'].map(df_conv.reset_index().set_index('systematic_name')['standard_name']).copy()
replacement_series = replacement_series.set_index('Gene')['standard_name'].to_dict()
df_hackett = df_hackett.replace({'Gene': replacement_series}).set_index('Gene')

df_hackett = df_conv.merge(df_hackett, how='left', left_index=True, right_index=True).copy()
df_hackett.set_index('systematic_name', inplace=True)

df_genes = df_genes.merge(df_hackett, how='left', left_index=True, right_index=True, sort=True)


###################################################
# KEMMEREN ET AL.
###################################################
df_kemmeren = pd.read_excel('../Data/Kemmeren2014.xlsx', engine='openpyxl', index_col=0, skiprows=[0,1], usecols=[1,3,5,6,8])
df_kemmeren.columns=['dFkh1_M', 'dFkh1_p-value', 'dFkh2_M', 'dFkh2_p-value']

def effect_dfkh(row):
    # Kemmeren et al. used a FC cutoff of 1.7 and p-value cutoff of 0.05
    # Use the above cutoff as a 'strong' up or down change
    # use a FC of 1.1 as medium threshold
    M = row[0]
    p = row[1]

    if p <= 0.05:
        verdict = 'Significant and '
    else:
        verdict = ''

    if abs(M) >= 100:
        verdict += 'suspicious value'
    elif M >= 0:
        if M >= np.log2(1.7):
            verdict += 'strongly up'
        elif M >= np.log2(1.1):
            verdict += 'up'
        else:
            verdict += 'weakly up'
    else:
        if M <= np.log2(1/1.7):
            verdict += 'strongly down'
        elif M <= np.log2(1/1.1):
            verdict += 'down'
        else:
            verdict += 'weakly down'

    return verdict.capitalize()

df_kemmeren['dFkh1'] = df_kemmeren[['dFkh1_M', 'dFkh1_p-value']].apply(effect_dfkh, axis=1)
df_kemmeren['dFkh2'] = df_kemmeren[['dFkh2_M', 'dFkh2_p-value']].apply(effect_dfkh, axis=1)

df_genes = df_genes.merge(df_kemmeren, how='left', left_index=True, right_index=True, sort=True)


###################################################
# VALIDATED?
###################################################
Fkh1_val = df_genes.query("not Fkh1_OE.isnull() and Fkh1_OE not in ('Invariant') \
                    and not dFkh1.isnull() \
                    and dFkh1 not in ('Weakly up','Weakly down','Suspicious value')", engine='python').index.tolist()
Fkh2_val = df_genes.query("not Fkh2_OE.isnull() and Fkh2_OE not in ('Invariant') \
                    and not dFkh2.isnull() \
                    and dFkh2 not in ('Weakly up','Weakly down','Suspicious value')", engine='python').index.tolist()
Fkh1_del = df_genes.query("not (not Fkh1_OE.isnull() and Fkh1_OE not in ('Invariant') ) \
                    and (not dFkh1.isnull() and dFkh1 not in ('Weakly up','Weakly down','Suspicious value') )", engine='python').index.tolist()
Fkh2_del = df_genes.query("not (not Fkh2_OE.isnull() and Fkh2_OE not in ('Invariant') ) \
                    and (not dFkh2.isnull() and dFkh2 not in ('Weakly up','Weakly down','Suspicious value') )", engine='python').index.tolist()
Fkh1_OE = df_genes.query("(not Fkh1_OE.isnull() and Fkh1_OE not in ('Invariant') ) \
                    and not (not dFkh1.isnull() and dFkh1 not in ('Weakly up','Weakly down','Suspicious value') )", engine='python').index.tolist()
Fkh2_OE = df_genes.query("(not Fkh2_OE.isnull() and Fkh2_OE not in ('Invariant') ) \
                    and not (not dFkh2.isnull() and dFkh2 not in ('Weakly up','Weakly down','Suspicious value') )", engine='python').index.tolist()

df_genes['Fkh1_validated'] = 'Neither'
df_genes['Fkh2_validated'] = 'Neither'
df_genes.loc[Fkh1_del, 'Fkh1_validated'] = 'Deletion only'
df_genes.loc[Fkh2_del, 'Fkh2_validated'] = 'Deletion only'
df_genes.loc[Fkh1_OE, 'Fkh1_validated'] = 'Overexpression only'
df_genes.loc[Fkh2_OE, 'Fkh2_validated'] = 'Overexpression only'
df_genes.loc[Fkh1_val, 'Fkh1_validated'] = 'Both'
df_genes.loc[Fkh2_val, 'Fkh2_validated'] = 'Both'


###################################################
# EXPORT EXCEL FILE
###################################################
# sort rows and columns
df_genes = df_genes[['standard_name', 'desc', 'expression_peak_phase', 'expression_peak_time', 'is_enzyme', 'KEGG_pathway', 
                    'Fkh1 MacIsaac', 'Fkh1 Venters', 'Fkh1 Ostrow', 'Fkh1 Mondeel', 'Fkh1 Rossi', 'Fkh1 Lupo z-score', 'Fkh1 Lupo',
                    'Fkh2 MacIsaac', 'Fkh2 Venters', 'Fkh2 Ostrow', 'Fkh2 Mondeel', 'Fkh2 Rossi', 'Fkh2 Lupo z-score', 'Fkh2 Lupo',
                    'Fkh1_times_target', 'Fkh2_times_target', 
                    'Fkh1_OE_0', 'Fkh1_OE_5', 'Fkh1_OE_8', 'Fkh1_OE_15', 'Fkh1_OE_30', 'Fkh1_OE_45', 'Fkh1_OE_60', 'Fkh1_OE_90', 'Fkh1_OE',
                    'Fkh2_OE_0', 'Fkh2_OE_5', 'Fkh2_OE_8', 'Fkh2_OE_15', 'Fkh2_OE_30', 'Fkh2_OE_45', 'Fkh2_OE_60', 'Fkh2_OE_90', 'Fkh2_OE',
                    'dFkh1_M', 'dFkh1_p-value', 'dFkh1', 'dFkh2_M', 'dFkh2_p-value', 'dFkh2',
                    'Fkh1_validated','Fkh2_validated']].sort_values(by='standard_name')

writer = pd.ExcelWriter("../Tables/Fkh1,2_regulation_summary.xlsx", engine='xlsxwriter')
workbook = writer.book

################### Full result sheet ###################
df_genes.to_excel(writer,sheet_name='Full table', index=True)


################### Fkh1 validated 4x+ targets ###################
df_4x_val = df_genes.query("Fkh1_times_target >= 4 and Fkh1_validated != 'Neither'")
df_4x_val.to_excel(writer,sheet_name='Fkh1_4x+_validated', index=True)


################### Fkh2 validated 4x+ targets ###################
df_4x_val = df_genes.query("Fkh2_times_target >= 4 and Fkh2_validated != 'Neither'")
df_4x_val.to_excel(writer,sheet_name='Fkh2_4x+_validated', index=True)


################### SET COLUMN WIDTH ###################
for worksheet_name in ['Full table','Fkh1_4x+_validated','Fkh2_4x+_validated']:
    worksheet = writer.sheets[worksheet_name]
    worksheet.set_column('A:B',12) # Systematic name - Standard name
    worksheet.set_column('C:C',50) # description
    worksheet.set_column('D:G',20) 
    worksheet.set_column('H:L',12) 
    worksheet.set_column('M:M',20) 
    worksheet.set_column('N:S',12) 
    worksheet.set_column('T:T',20) 
    worksheet.set_column('U:U',12) 
    worksheet.set_column('V:BB',15) 

    # freeze first row and column 1+2
    worksheet.freeze_panes(1, 2)

# SAVE
writer.save()