import numpy as np
import pandas as pd
import pdb
from os.path import join, dirname

### write rate matirx by HIV state
def write_rate_matrix_multiplier(df, risk_id):
    ### multiplier data 
    q_mat_blank_mlp = df.parse(sheet_name = 'counterfactual_q_mat{0}'.format(risk_id), index_col = 0).values 
    rates_indices_mlp = np.where(q_mat_blank_mlp != 0) 
    rate_data_mlp = df.parse(sheet_name = 'counterfactual_rate_data{0}'.format(risk_id), index_col = 0) # Number of HIV states x parameter number, col name: m1, m2,...
    rate_list_mlp = q_mat_blank_mlp[np.nonzero(q_mat_blank_mlp)]
    
    ### baseline data 
    q_mat_blank = df.parse(sheet_name = 'baseline_q_mat{0}'.format(risk_id), index_col = 0).values 
    rates_indices = np.where(q_mat_blank != 0) 
    rate_data = df.parse(sheet_name = 'baseline_rate_data{0}'.format(risk_id), index_col = 0)
    rate_list = q_mat_blank[np.nonzero(q_mat_blank)]

    # parameter
    num_age_group, num_disease_state = rate_data.shape[0], rate_data_mlp.shape[0]

    for state in range(num_disease_state):
        path = join(directory, 'rate_matrix_{0}_HIV{1}.xlsx'.format(risk_id, state))
        writer = pd.ExcelWriter(path, engine = 'xlsxwriter')
        for age_group in range(num_age_group):
            # calculate baseline transition rate matrix
            rate_matrix = np.zeros_like(q_mat_blank)
            rate_array = rate_data.loc[age_group, rate_list]
            rate_matrix[rates_indices] = rate_array
            
            # calculate counterfactual transition rate matrix based on disease status
            multiplier = np.ones_like(q_mat_blank_mlp)  
            rate_array_mlp = rate_data_mlp.loc[state, rate_list_mlp]
            multiplier[rates_indices_mlp] = rate_array_mlp
            rate_matrix_mlp = np.multiply(rate_matrix, multiplier)
            
            ### pay attention to here 
            ### it subjects to change
            # indicate it is equal to multiplier itself
            # rate_matrix_mlp[-4:-1,-1] = multiplier[-4:-1,-1] # modify cancer state mortality: local, regional, distant
            df_new = pd.DataFrame (rate_matrix_mlp)
            df_new.to_excel(writer, sheet_name = "age_group{0}".format(age_group))
        try:
            writer.close()
        except Exception as e:
            print(f"Error occurred: {e}")

if __name__ == "__main__":
    ## remember to pay attention to the part where mortality was modified
    current_dir = dirname(__file__)
    directory = './rate_matrix_baseline'
    df = pd.ExcelFile(join(current_dir,"./rate_matrix_baseline.xlsx"))
    print("current working directory:", directory)
    num_risk_group = 3
    for i in range(num_risk_group):
        write_rate_matrix_multiplier(df, i)



   
    
    