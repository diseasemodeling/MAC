import numpy as np
import csv
import pdb
import pickle

global degree_dist_bin_all
degree_dist_bin_all = np.asarray([0, 1, 2, 4, 8, 16, 32, 64, 128])

global age_lb
age_lb = np.asarray([13,17,24,29,34,39,44,65])
global age_ub
age_ub = np.asarray([17,24,29,34,39,44,65,150])


def find_age_index(age):
    i = 0
    found = False
    age_index = -1
    while i < len(age_ub) and found == False:
        if age <= age_ub[i]:
            age_index = i
            found = True
        i += 1
    ### just in case 
    if age > age_ub[-1]: 
        age_index = -1
        # count_error_age[index] += 1
    elif age < age_lb[0]: 
        age_index = 0
        # count_error_age[index] += 1
    return age_index

### checked 
def find_degree_bin(degree):
    column_bin = -1
    if degree <= degree_dist_bin_all[-1] and degree >= degree_dist_bin_all[0]:
        column_bin = np.where(degree_dist_bin_all >= degree)[0][0]
    # else: 
    #     raise ValueError ("Degree is not within the range")
    return column_bin

def return_actual_model_matrix(age_l, degree_l, sex_l, d1_l, d2_l, model_matrix, risk_ind):
    actual_model_matrix = np.zeros_like(model_matrix[risk_ind])
    for m in range(len(degree_l)):
        k = int(sex_l[m])
        if k == risk_ind:
            i = find_age_index(float(age_l[m]))
            j = find_degree_bin(float(degree_l[m]))

            p = int(d1_l[m])
            q = int(d2_l[m])
            actual_model_matrix[p, i, j, 0, q]+= 1
    actual_model_matrix[0]=model_matrix[risk_ind][0]
    return actual_model_matrix


def read_file(file_name):
    with open(file_name, "rb") as f:
        debug = pickle.load(f)
        model_matrix = debug['model']
        model_change = debug['model_change']
        age_l = np.array(debug['age'])
        degree_l = np.array(debug['deg'])
        sex_l = np.array(debug['sex'])
        d1_l = np.array(debug['d1'])
        d2_l = np.array(debug['d2'])
    return model_matrix, model_change,age_l, degree_l, sex_l, d1_l, d2_l


def check_model_matrix_consistency(age_l, degree_l, sex_l, d1_l, d2_l, model_matrix):
    """
    Check the consistency of the model matrix with the agent data.
    """
    # Convert lists to numpy arrays
    age_l = np.array(age_l)
    degree_l = np.array(degree_l)
    sex_l = np.array(sex_l)
    d1_l = np.array(d1_l)
    d2_l = np.array(d2_l)
    
    # Initialize matrices
    num_risk_group = 3
    actual_model_matrix = {risk_ind: np.zeros_like(model_matrix[risk_ind]) for risk_ind in range(num_risk_group)}
    differences = {risk_ind: np.zeros_like(model_matrix[risk_ind]) for risk_ind in range(num_risk_group)}

    # Calculate the actual model matrix
    for risk_ind in range(num_risk_group):
        actual_model_matrix[risk_ind] = return_actual_model_matrix(age_l, degree_l, sex_l, d1_l, d2_l, model_matrix,risk_ind)
        differences[risk_ind] = actual_model_matrix[risk_ind] - model_matrix[risk_ind]

    # Check for inconsistencies
    inconsistencies = {risk_ind: len(np.where(differences[risk_ind] != 0)[0]) > 0 for risk_ind in range(num_risk_group)}
    inconsistencies_list = [inconsistencies[risk_ind] for risk_ind in range(num_risk_group)]
    return inconsistencies_list

        


# model_matrix, model_change, age_l, degree_l, sex_l, d1_l, d2_l = read_file("debug109.pkl")

# actual_model_matrix = {}
# num_risk_group = 3
# for risk_ind in range(num_risk_group):
#     actual_model_matrix[risk_ind] = return_actual_model_matrix(age_l, degree_l, sex_l, d1_l, d2_l, model_matrix,risk_ind)
#     dif = actual_model_matrix[risk_ind]- model_matrix[risk_ind]
#     print(risk_ind)
#     if np.any(dif>0): print('dif>0', np.where(dif>0))
#     if np.any(dif<0): print('dif<0', np.where(dif<0))    



# pdb.set_trace()