import pandas as pd
import numpy as np
def convert():
    raw_data = pd.read_csv('/Users/linze/PycharmProjects/untitled/datacleaning/selected_labeled.csv')
    raw_data = raw_data.transpose ()
    gene_name = raw_data.iloc[0,]
    raw_data = raw_data.drop (["Unnamed: 0"])
    raw_data = raw_data.rename (gene_name, axis='columns')
    dmm_1 = raw_data.loc[raw_data['Label'] == 0]
    dmm_1 = dmm_1.drop (columns=['Label'])
    sham_1 = raw_data.loc[raw_data['Label'] == 1]
    sham_1 = sham_1.drop (columns=['Label'])
    dmm_6 = raw_data.loc[raw_data['Label'] == 2]
    dmm_6 = dmm_6.drop (columns=['Label'])
    sham_6 = raw_data.loc[raw_data['Label'] == 3]
    sham_6 = sham_6.drop (columns=['Label'])
    gene_name=gene_name.drop(index=1291)
    dmm_log = np.empty([16,dmm_1.shape[1]])
    index=0
    for index_1, row_1 in dmm_1.iterrows ():
        for index_6, row_6 in dmm_6.iterrows ():
            tmp_1=row_1.to_numpy(dtype=float)
            tmp_6=row_6.to_numpy(dtype=float)
            tem_result=np.divide(tmp_1,tmp_6)
            dmm_log[index]=tem_result
            index+=1
    dmm_log = pd.DataFrame (data=dmm_log,columns=gene_name)
    sham_log = np.empty([16,dmm_1.shape[1]])
    index1=0
    for index_1, row_1 in sham_1.iterrows ():
        for index_6, row_6 in sham_6.iterrows ():
            tmp_1=row_1.to_numpy(dtype=float)
            tmp_6=row_6.to_numpy(dtype=float)
            tem_result=np.log2(np.divide(tmp_1,tmp_6))
            sham_log[index1]=tem_result
            index1+=1
    sham_log = pd.DataFrame (data=sham_log,columns=gene_name)
    result=pd.concat([dmm_log,sham_log])
    label=['dmm', 'dmm', 'dmm', 'dmm', 'dmm', 'dmm', 'dmm', 'dmm', 'dmm', 'dmm', 'dmm', 'dmm', 'dmm', 'dmm', 'dmm', 'dmm', 'sham', 'sham', 'sham', 'sham', 'sham', 'sham', 'sham', 'sham', 'sham', 'sham', 'sham', 'sham', 'sham', 'sham', 'sham', 'sham']
    result['label']=label
    print(result.head())
    result.to_csv("log2_ratio_modify.csv")



# for each_6 in dmm_6:


if __name__ == "__main__":
    convert()
