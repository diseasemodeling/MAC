import numpy as np
from numba import njit,prange
import pandas as pd
import pickle
import pdb
np.set_printoptions(suppress=True)
from .config_reader import ConfigReader
from .sexual_behavior import SexualBehavior
from .utils import config_path, data_dir

config = ConfigReader(config_path)


TRANS_PROB_KEY = "trans_prob"
TRANSMISSION_PAIR_KEY = "transmission_pair_index"


DISEASE_FOLDER_NAME_KEY = "disease_folder_name"
NUM_DISEASE_STATE_KEY = "num_disease_state"
BASE_DISEASE_KEY = "base_disease_id"
MODEL_DISEASE_KEY = "modeling_disease_id"
INCREASED_RISK_KEY = "increased_risk"
SUBGROUP_KEY = "num_subgroup_type"

INFECTION_AREA_INDEX_KEY = "infection_area_index"
INFECTION_TYPE_INDEX_KEY = "infection_type_index"
DISEASE_STATE_INDEX_KEY = "disease_state_index"
DIAGNOSED_STATE_INDEX_KEY = "diagnosed_state_index"
DUMMY_STATE_INDEX_KEY = "dummy_state_index"
UNKNOWN_INDEX_KEY = "unknown_index"
NUM_DIAGNOSED_STATE_KEY = "num_diagnosed_state"

SUSCEPTIBLE_INDEX_KEY = "susceptible_index"
INFECTED_INDEX_KEY = "infected_index"
IMMUNITY_INDEX_KEY = "immunity_index"
PRECANCER_INDEX_KEY = "precancer_index"
CANCER_INDEX_KEY = "cancer_index"
VACCINE_INDEX_KEY = "vaccine_index"
DEATH_INDEX_KEY = "death_state_index"
NONAGENTS_INDEX_KEY = "nonagents_index"

CANCER_MORTALITY_RATE_KEY = "cancer_mortality_rate"
NATURAL_MORTALITY_RATE_KEY = "natural_mortality_rate"
INIT_PROP_INF_KEY = "init_prop_inf"

SCREEN_PRECANCER_KEY = "screening_rate_precancer"
SCREEN_CANCER_KEY = "screening_rate_cancer"
SCREEN_COMPLIANCE_KEY = "screening_compliance"
SCREEN_SENSITIVITY_KEY = "screening_sensitivity"
FOLLOWUP_SCREEN_KEY = "follow_up_screening_compliance"
FOLLOWUP_TREATMENT_KEY = "follow_up_treatment_compliance"
FOLLOWUP_SCREEN_SENSITIVITY_KEY = "follow_up_screening_sensitivity"
PROP_TREATMENT_KEY = "prop_treatment"
VACCINE_EFFICACY_KEY = "vaccine_efficacy"
VACCINE_SHEET_KEY =  "vaccine_rate_name"

def return_index_list(df, col_num, state_num):
    index = df[df. iloc[:, col_num] == state_num].index.tolist()
    return [int(i) for i in index]  # list

def parse_list2numpy(list):
    return np.array(list)


# @njit(parallel=True)
# def disease_progression_multinomial(model_change, model_matrix, transition_prob_matrix, num_disease_state_base, num_age_group, num_degree_bin, num_jurisdictions, num_disease_state_modeling, master_seed):
#     for i in prange(num_disease_state_base):       
#     # for i in prange(num_disease_state_base):
#         for j in prange(num_age_group):
#             for k in prange(num_degree_bin):
#                 for m in prange(num_jurisdictions):
#                     for n in prange(num_disease_state_modeling):
#                         n_val = np.int64(np.floor(model_matrix[i,j,k,m,n]))
#                         for p in range(num_disease_state_modeling):
#                             if n!=p:
#                                 pvals = transition_prob_matrix[i,j,k,m,n,p]
#                                 sample = np.random.binomial(n_val, pvals)
#                                 model_change[i,j,k,m,n,p] = sample
#                                 n_val = max(n_val - sample, 0)

def disease_progression_multinomial(model_change, model_matrix, transition_prob_matrix, num_disease_state_base, num_age_group, num_degree_bin, num_jurisdictions, num_disease_state_modeling, master_seed):
    np.random.seed(master_seed)
    for i in range(1, num_disease_state_base):       # 1 refers to self.nonagents_state 
        for j in range(num_age_group):
            for k in range(num_degree_bin):
                for m in range(num_jurisdictions):
                    for n in range(num_disease_state_modeling):
                        n_val = np.int64(np.floor(model_matrix[i,j,k,m,n]))
                        pvals = transition_prob_matrix[i,j,k,m,n]
                        model_change[i,j,k,m,n] = np.random.multinomial(n_val, pvals)
                        
                                        
class Disease():
    def __init__(self, 
                 sb = None,
                 model_matrix = {},
                 model_change = {},
                 netlogo_change = {},
                 nonagents_state = 0,
                 time_unit = 2, 
                 mode = None,
                 increased_risk = 0):
        self.mode = mode
        self.time_unit = time_unit
        self.increased_risk = increased_risk
        
        self.model_matrix = model_matrix
        self.model_change = model_change
        self.netlogo_change = netlogo_change
        self.nonagents_state = nonagents_state

        
        if sb is None:
            self.sb = SexualBehavior(time_unit = time_unit)
        else:
            self.sb = sb
        self.num_age_group =  self.sb.num_age_group
        self.num_degree_bin = self.sb.num_degree_bin
        self.num_jurisdictions = self.sb.num_jurisdictions
        self.num_risk_group = self.sb.num_risk_group
        self.num_adjusted_ptnr = self.sb.num_adjusted_ptnr
        
        self.trans_prob = parse_list2numpy(config.get(TRANS_PROB_KEY))[:self.num_risk_group] 
        self.calc_transmission_ptnrship_probability()
        
        self.num_disease_state = config.get(NUM_DISEASE_STATE_KEY)
        self.base_disease_id = config.get(BASE_DISEASE_KEY)
        self.modeling_disease_id = config.get(MODEL_DISEASE_KEY)
        self.num_subgroup_type = config.get(SUBGROUP_KEY)
        self.death_state = config.get(DEATH_INDEX_KEY)
        self.init_prop_inf = config.get(INIT_PROP_INF_KEY)
        
        disease_folder_name = config.get(DISEASE_FOLDER_NAME_KEY)
        
        self.initialize_variables()
        self.read_intervention_data()
        self.read_natural_hist_data(disease_folder_name)
        self.calc_base_transition_rate_matrix(disease_folder_name)
        self.define_disease_state_index()
        

    @staticmethod
    def get_num_states():
        num_disease_state = config.get(NUM_DISEASE_STATE_KEY)
        base_disease_id = config.get(BASE_DISEASE_KEY)
        modeling_disease_id = config.get(MODEL_DISEASE_KEY)
        return num_disease_state[base_disease_id], num_disease_state[modeling_disease_id]

    def initialize_variables(self):
        self.transition_rate_matrix = {}
        self.transition_prob_matrix = {}
        self.prop_inf_by_group = np.zeros((self.num_risk_group, self.num_age_group,\
                                        self.num_degree_bin, self.num_jurisdictions, \
                                         self.num_subgroup_type))
    # read natural history data and multiplier data 
    ## checked ##
    def read_natural_hist_data(self, disease_folder_name):
        rate_matrix_name = disease_folder_name + ".xlsx"
        
        # print("current rate data directory: ", disease_folder_name)
        # print("current data file: ", rate_matrix_name)
        
        excel_path = data_dir/ rate_matrix_name
        df = pd.ExcelFile(excel_path)
        self.rate_data = {n: df.parse(sheet_name = 'baseline_rate_data{0}'.format(n), index_col = 0) for n in range(self.num_risk_group)}
        self.multiplier_data = {n: df.parse(sheet_name = 'counterfactual_rate_data{0}'.format(n), index_col = 0) for n in range(self.num_risk_group)} 
        self.base_transition_rate_matrix = {n: np.zeros((self.num_disease_state[self.base_disease_id][n], self.num_age_group,\
                                            self.num_disease_state[self.modeling_disease_id][n], self.num_disease_state[self.modeling_disease_id][n]))for n in range(self.num_risk_group)}
        self.q_mat = {n: df.parse(sheet_name = 'baseline_q_mat{0}'.format(n), index_col = 0).values for n in range(self.num_risk_group)}
        # state ind 
        self.state_list = {n: df.parse(sheet_name = 'state{0}'.format(n), index_col = 0) for n in range(self.num_risk_group)}
    
    # read intervention related data
    def read_intervention_data(self):
        excel_path = data_dir/ config.get("intervention_data")
        df = pd.ExcelFile(excel_path)
        self.vaccine_rate = df.parse(sheet_name = config.get(VACCINE_SHEET_KEY), index_col = 0).values 
        self.screening_rate_precancer =df.parse(sheet_name =config.get(SCREEN_PRECANCER_KEY)).values
        self.screening_rate_cancer = df.parse(sheet_name =config.get(SCREEN_CANCER_KEY)).values
        self.screening_sensitivity = df.parse(sheet_name =config.get(SCREEN_SENSITIVITY_KEY)).values
        self.follow_up_screening_compliance = df.parse(sheet_name =config.get(FOLLOWUP_SCREEN_KEY)).values
        self.prop_treatment = df.parse(sheet_name =config.get(PROP_TREATMENT_KEY)).values
        
        self.vaccine_efficacy = config.get(VACCINE_EFFICACY_KEY)
        self.screening_compliance = config.get(SCREEN_COMPLIANCE_KEY)
        self.follow_up_treatment_compliance = config.get(FOLLOWUP_TREATMENT_KEY)
        self.follow_up_screening_sensitivity = config.get(FOLLOWUP_SCREEN_SENSITIVITY_KEY)
    
    ## calculate basic transition rate matrix 
    def calc_base_transition_rate_matrix(self, disease_folder_name):    
        for risk_ind in range(self.num_risk_group):
            for CD4_state in range(self.num_disease_state[self.base_disease_id][risk_ind]):
                ### modify increased risk 
                if self.increased_risk == 0 or self.increased_risk == 3: CD4_state_mod = 0 
                else: CD4_state_mod = CD4_state
                excel_path = data_dir/ disease_folder_name/'rate_matrix_{0}_HIV{1}.xlsx'.format(risk_ind, CD4_state_mod)
                data = pd.ExcelFile(excel_path)
                for age_group in range(self.num_age_group):
                    # if CD4_state>=1:
                    #     self.base_transition_rate_matrix[risk_ind][CD4_state, age_group] = data.parse(sheet_name = "age_group%s"%age_group, header = 0, index_col = 0).values *2
                    # else:
                    self.base_transition_rate_matrix[risk_ind][CD4_state, age_group] = data.parse(sheet_name = "age_group%s"%age_group, header = 0, index_col = 0).values
            
            # modify natural mortality for cancer stage   
            self.modify_mortality_parameters(risk_ind) 

            
    # broadcast transition rate matrix to other index
    def transform_transition_rate_matrix(self, risk_ind):
        expanded_matrix1 = self.base_transition_rate_matrix[risk_ind][:, :, np.newaxis, :, :]
        reshaped_matrix1 = np.repeat(expanded_matrix1, self.num_degree_bin, axis=2)
        expanded_matrix2 =  reshaped_matrix1 [:, :, :, np.newaxis, :, :]
        final_matrix = np.repeat(expanded_matrix2, self.num_jurisdictions, axis=3)
        self.transition_rate_matrix[risk_ind] = final_matrix   
   
    
    def update_screening_parameters(self, risk_ind, screening_condition = False):
        if (risk_ind == 0 or risk_ind == 2) and screening_condition:
            num_precancer_state = len(self.precancer_index[risk_ind])
            
            # Compute base rates for screening precancer and cancer
            base_precancer = self.screening_rate_precancer * self.screening_compliance * self.follow_up_screening_sensitivity * self.follow_up_treatment_compliance
            base_cancer = self.screening_rate_cancer * self.screening_compliance * self.follow_up_screening_sensitivity
            
            # Compute the rate of screening for precancer and cancer states using vectorized operations
            r_screening_precancer = np.outer(base_precancer, self.screening_sensitivity[:,:num_precancer_state]) * self.follow_up_screening_compliance[:,:num_precancer_state] * self.prop_treatment
            r_screening_cancer = np.outer(base_cancer, self.screening_sensitivity[:,num_precancer_state:]) * self.follow_up_screening_compliance[:,num_precancer_state:]
            r_screening = np.hstack((r_screening_precancer, r_screening_cancer))
            
            start_id = 6 # hard code for initial disease state CIN 1 or AIN 1

            for id in range(start_id, start_id+r_screening.shape[1]):
                dif = start_id + r_screening_precancer.shape[1] - id 
                if dif > 0:
                    state = [i for i in return_index_list(self.state_list[risk_ind], 2, id)]
                    transition_id = [i[1] for i in list(self.infection_index[risk_ind].values())]
                else:
                    state = self.cancer_index[risk_ind][0][dif]
                    transition_id = self.cancer_index[risk_ind][1][dif]
                
                temp_var = self.transition_rate_matrix[risk_ind][...,state,:]
                next_temp_var = temp_var[...,transition_id]
                mask = next_temp_var > 0 
                if dif > 0:
                    created_arr = np.repeat(np.repeat(np.repeat(np.repeat(np.repeat(r_screening[:,id-start_id][None,...], next_temp_var.shape[0], axis =0)[...,None], next_temp_var.shape[2],axis = 2)[...,None], next_temp_var.shape[3],axis =3)[...,None], next_temp_var.shape[4],axis =4)[...,None],len(transition_id),axis =5)
                else:
                    created_arr = np.repeat(np.repeat(np.repeat(r_screening[:,id-start_id][None,...], next_temp_var.shape[0], axis =0)[...,None], next_temp_var.shape[2],axis = 2)[...,None], next_temp_var.shape[3],axis =3)
                next_temp_var[mask] += created_arr[mask]
                temp_var[...,transition_id] = next_temp_var
                self.transition_rate_matrix[risk_ind][...,state,:] = temp_var

    
    # Add natural mortality rate on top of cancer related mortality
    # checked
    def modify_mortality_parameters(self, risk_ind):
        cancer_mortality_rate = config.get(CANCER_MORTALITY_RATE_KEY)
        natural_mortality_rate = config.get(NATURAL_MORTALITY_RATE_KEY)
       
        for rate_ind in cancer_mortality_rate:
            rate_values = np.array(self.rate_data[risk_ind][natural_mortality_rate].values)[None, :, None]
            condition = self.q_mat[risk_ind] == rate_ind
            self.base_transition_rate_matrix[risk_ind][..., condition] += rate_values
                      
    ## calculate transmission probability per partnership                
    def calc_transmission_ptnrship_probability(self):
        num_ptnrs = np.copy(self.sb.num_ptnrs)
        # num_ptnrs = np.copy(self.sb.num_ptnrs/self.time_unit)
        num_acts = np.copy(self.sb.num_acts/self.time_unit)
        
        mask_main_ptnr = self.sb.num_ptnrs_casual  <= 1    # determine if main/ casual partnerships by monthly number of partners 
        # mask_main_ptnr = self.sb.num_ptnrs/12  <= 1 
        prop_protected_acts = mask_main_ptnr * self.sb.prop_protected_acts_main[:,:,None] \
                             + ~mask_main_ptnr * self.sb.prop_protected_acts_casual[:,:,None]

        np.seterr(divide='ignore', invalid='ignore')
        num_exposure = np.where(num_ptnrs == 0, 0,
                                np.where(
                                    num_ptnrs >= 1, 
                                    num_acts[:, :, None] / num_ptnrs, 
                                    num_acts[:, :, None]
                                )
                                ) 
        # num_exposure = np.where(num_ptnrs > 0, num_acts[:,:,None]/num_ptnrs, 0)
        not_trans_protected = (1- ((1-self.sb.condom_efficacy[:, None, None]) * \
                                 self.trans_prob[:, None, None, :]))**(prop_protected_acts[:,:, :, None]\
                                  *  num_exposure[:,:, :, None])
        not_trans_unprotected = (1- self.trans_prob[:,None, None,:])**\
                                  ((1-prop_protected_acts[:,:,:, None]) * num_exposure[:,:,:,None])
        self.transmission_prob_ptnrship = 1 - not_trans_protected * not_trans_unprotected
    
    
    def transmission(self):
        self.calc_prevalence()
        self.calc_FOI()
    
    def transition(self, 
                   screening_condition = False, 
                   vaccination_condition = False, 
                   year = 2006,
                   random_seed = 1076):
        for risk_ind in range(self.num_risk_group):
            
            # broadcast transition rate matrix to degree bin and jurisdiction
            self.transform_transition_rate_matrix(risk_ind)
            
            self.update_transition_rate_matrix(risk_ind, screening_condition, vaccination_condition, year)
            
            self.update_transition_prob_matrix(risk_ind)
           
            self.wrapper_function(risk_ind, random_seed)
           
            self.update_transition(risk_ind)

        
    def update_transition(self, risk_ind):
        if self.mode == "netlogo":
            self.model_change[risk_ind][self.nonagents_state] =  self.model_matrix[risk_ind][self.nonagents_state,...,None] * self.transition_rate_matrix[risk_ind][self.nonagents_state]
            self.model_matrix[risk_ind][self.nonagents_state] += self.model_change[risk_ind][self.nonagents_state].sum(-2)
        else:
            self.model_change[risk_ind][self.nonagents_state] =  self.model_matrix[risk_ind][self.nonagents_state,...,None] * self.transition_rate_matrix[risk_ind][self.nonagents_state]
            self.model_change[risk_ind][1:-1] =  np.floor(self.model_matrix[risk_ind][1:-1,...,None] * self.transition_rate_matrix[risk_ind][1:-1]) 
            self.model_matrix[risk_ind] += self.model_change[risk_ind].sum(-2)

    ### calculate prevalence by type
    ## checked 
    def calc_prevalence(self):
        for risk_ind in range(self.num_risk_group):
            death_state_base = self.death_state[self.base_disease_id]
            death_state_model = self.death_state[self.modeling_disease_id][risk_ind]

            # Calculate population size for the current group
            total_popsize_by_group = np.delete(self.model_matrix[risk_ind],death_state_base,axis =0)
            total_popsize_by_group = np.delete(total_popsize_by_group,death_state_model,axis =-1)
            popsize_by_group = total_popsize_by_group.sum((0,-1))

            for type_ind in range(self.num_subgroup_type):
                # Calculate the number of infections for the current group and type
                type = self.infection_indexBYTYPE[risk_ind][type_ind]
                num_inf_by_group = total_popsize_by_group[...,type].sum((0,-1))
 
                # Handle divide by zero or invalid values
                np.seterr(divide='ignore', invalid='ignore')

                # Calculate and store prevalence for the current group and type
                condition = popsize_by_group >= 1
                prevalence = num_inf_by_group / popsize_by_group
                self.prop_inf_by_group[risk_ind, ..., type_ind] = np.where(condition, prevalence, 0)

    
    ### calculate force of infection 
    def calc_FOI(self):
        # Calculate the number of effective partners
        num_effective_ptnr = (
            self.sb.num_adjusted_ptnr[..., None] *
            self.prop_inf_by_group[None, :, None, :, None, :, None, :, :]
        )

        # Calculate the effective transmission rate
        effective_trans_rate = (
            num_effective_ptnr *
            self.transmission_prob_ptnrship[:, None, :, None, :, None, None, None, :]
        )

        # Sum along specified axes to compute the force of infection
        self.force_of_infection = effective_trans_rate.sum((1, 3, 5, 7))
    
    ## update transition rate matrix
    def update_transition_rate_matrix(self, risk_ind, screening_condition, vaccination_condition, year):
        self.update_screening_parameters(risk_ind,screening_condition)
        self.update_vaccination_parameters(risk_ind, vaccination_condition, year) 
        self.transition_rate_matrix[risk_ind] *= 1/self.time_unit
        self.update_transmission_parameters(risk_ind)
              
        row_sum = self.transition_rate_matrix[risk_ind].sum(-1)
        idx = np.arange(self.num_disease_state[self.modeling_disease_id][risk_ind])
        self.transition_rate_matrix[risk_ind][...,idx,idx] = -1 * row_sum 

    
    # update transition probability matrix 
    def update_transition_prob_matrix(self, risk_ind):
        idx = np.arange(self.num_disease_state[self.modeling_disease_id][risk_ind])
        # sum of the rates among other events that will lead to departure of the compartment 
        row_sum = -self.transition_rate_matrix[risk_ind][...,idx,idx]
        # division of the rates that will lead to departure of the compartment among sum of the rates among all events
        div = self.transition_rate_matrix[risk_ind]/row_sum[...,None] 
        exp = 1-np.exp(-row_sum[...,None])
        self.transition_prob_matrix[risk_ind] = div * exp
        self.transition_prob_matrix[risk_ind][...,idx,idx] = np.exp(-row_sum)
        self.transition_prob_matrix[risk_ind] = np.nan_to_num(self.transition_prob_matrix[risk_ind])
        ### manually normalization to 1
        self.transition_prob_matrix[risk_ind] /=self.transition_prob_matrix[risk_ind].sum(-1)[...,None]
        

    ## update force of infection 
    def update_transmission_parameters(self, risk_ind):
        # Determine whether transmission rate would be increased
        no_trans_multiplier_condition = self.increased_risk % 2 == 0
        transmission_pair_index = config.get(TRANSMISSION_PAIR_KEY)

        # Add one dimension for base disease 
        expanded_matrix = self.force_of_infection[:, None, ...]
        lambda_v = np.repeat(expanded_matrix, self.num_disease_state[self.base_disease_id][risk_ind], axis=1)

        for i, pair in enumerate(transmission_pair_index):
            # Define the multiplier
            val = lambda_v[risk_ind, ..., i] 
            mul = np.array(self.multiplier_data[risk_ind][pair[1]])
            expanded_mul = mul[:, None, None, None]
            multiplier = val if no_trans_multiplier_condition else val * expanded_mul

            # Update the transition_rate_matrix based on the pair index
            self.transition_rate_matrix[risk_ind][..., self.q_mat[risk_ind] == pair[0]] = multiplier[...,None]    
    
    def wrapper_function(self, risk_ind, random_seed):
        num_disease_state_base = self.num_disease_state[self.base_disease_id][risk_ind]
        num_disease_state_modeling = self.num_disease_state[self.modeling_disease_id][risk_ind]
        num_age_group = self.num_age_group
        num_degree_bin = self.num_degree_bin
        num_jurisdictions = self.num_jurisdictions
        
        self.model_change[risk_ind] = np.zeros_like(self.transition_prob_matrix[risk_ind])
        # self.netlogo_change[risk_ind] = np.zeros_like(self.transition_prob_matrix[risk_ind])

        disease_progression_multinomial(self.model_change[risk_ind], 
                                        self.model_matrix[risk_ind], 
                                        self.transition_prob_matrix[risk_ind], 
                                        num_disease_state_base, 
                                        num_age_group, 
                                        num_degree_bin, 
                                        num_jurisdictions, 
                                        num_disease_state_modeling, 
                                        random_seed)
        
        # mannually calculate change in model matrix for the persons not leaving the compartment 
        idx = np.arange(self.num_disease_state[self.modeling_disease_id][risk_ind])
        self.model_change[risk_ind][...,idx,idx] = 0 
        row_sum = self.model_change[risk_ind].sum(-1) 
        self.model_change[risk_ind][...,idx,idx] = np.minimum(-1 * row_sum,  self.model_matrix[risk_ind][...,idx])
            
    ### to check###
    def update_vaccination_parameters(self, risk_ind, vaccination_condition, year):        
        if vaccination_condition: 
            initial_year = 2006
            r_vaccination = self.vaccine_rate[year-initial_year,risk_ind] * self.vaccine_efficacy  
            first_age_group_index = 0
            if year < 2015:
                self.transition_rate_matrix[risk_ind][:,first_age_group_index,...,self.susceptible_index[risk_ind], self.vaccination_index[risk_ind][0]] = r_vaccination
            else: 
                self.transition_rate_matrix[risk_ind][:,first_age_group_index,...,self.susceptible_index[risk_ind], self.vaccination_index[risk_ind][1]] = r_vaccination
        
    # initialize infection to population
    ## assume equally distributed by subgroup type
    def initialize_infection2population(self):
        for risk_ind in range(self.num_risk_group):
            # model_matrix_slice = self.model_matrix[risk_ind, self.nonagents_state,:,self.non_zero_degree_bin:]
            model_matrix_slice = self.model_matrix[risk_ind][self.nonagents_state]
            pop_to_change = model_matrix_slice[...,self.susceptible_index[risk_ind]] * self.init_prop_inf
            model_matrix_slice[...,self.susceptible_index[risk_ind]] -= pop_to_change
            
            temp_inf_ind = np.array(list(self.infection_index[risk_ind].values()))[:, 0]
            model_matrix_slice[..., temp_inf_ind[:self.num_subgroup_type]] += pop_to_change/self.num_subgroup_type
    
    #### checked
    def define_disease_state_index(self):
        susceptible_index = config.get(SUSCEPTIBLE_INDEX_KEY)
        infected_index = config.get(INFECTED_INDEX_KEY)
        immunity_index = config.get(IMMUNITY_INDEX_KEY)
        precancer_index = config.get(PRECANCER_INDEX_KEY)
        cancer_index = config.get(CANCER_INDEX_KEY)
        vaccine_index = config.get(VACCINE_INDEX_KEY) 

        infection_area_index = config.get(INFECTION_AREA_INDEX_KEY)
        infection_type_index = config.get(INFECTION_TYPE_INDEX_KEY)
        disease_index = config.get(DISEASE_STATE_INDEX_KEY)
        diagnosed_index = config.get(DIAGNOSED_STATE_INDEX_KEY)
        dummy_index = config.get(DUMMY_STATE_INDEX_KEY)
        unknown_index = config.get(UNKNOWN_INDEX_KEY)
        num_diagnosed_state = config.get(NUM_DIAGNOSED_STATE_KEY)
        # self.death_state = [i - 1 for i in self.num_disease_state]
        
        # vaccination index
        self.vaccination_index = {n: vaccine_index for n in range(self.num_risk_group)}
        
        ## overall hpv susceptible
        self.susceptible_index = {n: [i for k in susceptible_index 
                                      for i in return_index_list(self.state_list[n], disease_index, k)]  
                                      for n in range(self.num_risk_group)}
        
        # natural immunity 
        self.natural_immunity_index =  {n: [i for k in immunity_index 
                                            for i in return_index_list(self.state_list[n], disease_index, k)] 
                                            for n in range(self.num_risk_group)}
        
        # HPV infection 
        self.infection_index = {n: {j: [i for k in infected_index 
                                        for i in return_index_list(self.state_list[n], disease_index, k) 
                                        if i in return_index_list(self.state_list[n], infection_type_index , j+1)]  
                                        for j in range(self.num_subgroup_type)} for n in range(self.num_risk_group)}
        # overall hpv infection 
        self.overall_infection_index = {n: [i for i in range(self.num_disease_state[self.modeling_disease_id][n]) 
                                             if i in return_index_list(self.state_list[n], diagnosed_index, unknown_index) 
                                             if i not in self.vaccination_index[n]
                                             if i not in self.natural_immunity_index[n]
                                             if i not in self.susceptible_index[n]]
                                             for n in range(self.num_risk_group)}
        # cancer_index 
        self.cancer_index = {n: {j: [i for i in return_index_list(self.state_list[n], diagnosed_index, j+1)] 
                                       for j in range(num_diagnosed_state)}
                                       for n in range(self.num_risk_group)}
        # precancer_index 
        self.precancer_index = {n: {j: [i for k in precancer_index 
                                          for i in return_index_list(self.state_list[n], disease_index, k) 
                                          if i in return_index_list(self.state_list[n], infection_type_index, j+1)]  
                                          for j in range(self.num_subgroup_type) } 
                                          for n in range(self.num_risk_group)}
        
        # overall hpv infection by type
        self.infection_indexBYTYPE = {n: {j: [i for i in return_index_list(self.state_list[n], infection_type_index, (j+1)) 
                                                if i not in self.natural_immunity_index[n]]
                                                for j in range(self.num_subgroup_type)} 
                                                for n in range(self.num_risk_group)}
        
    # calculate genotype frequency 
    def calc_genotype_frequency(self, risk_ind):
        num_precancer_state = len(self.precancer_index[risk_ind])
        # Initialize arrays for genotype frequencies
        genotype_frequency_CIN = np.zeros((num_precancer_state, self.num_subgroup_type, self.num_age_group))
        genotype_frequency_normal = np.zeros((self.num_subgroup_type, self.num_age_group))
       
        model_mat_risk = np.delete(self.model_matrix[risk_ind],self.death_state[self.base_disease_id],axis =0)

        for i in range(self.num_subgroup_type):
            # Genotype frequency among CIN 1-3
            genotype_frequency_CIN[:, i] = model_mat_risk[..., self.precancer_index[risk_ind][i]].sum((0, 2, 3)).T
            genotype_frequency_normal[i] = model_mat_risk[..., self.infection_index[risk_ind][i]].sum((0, 2, 3, 4))

        # Genotype frequency among all population
        genotype_frequency_general = (genotype_frequency_normal + genotype_frequency_CIN.sum(0))

        # Normalize frequencies and convert to percentages
        total_for_normalization_CIN = genotype_frequency_CIN.sum(1)[:, None, :]
        genotype_frequency_CIN = np.where(total_for_normalization_CIN != 0,(genotype_frequency_CIN / total_for_normalization_CIN) * 100,0)
        
        total_for_normalization_normal = (model_mat_risk[..., self.susceptible_index[risk_ind]].sum((0, 2, 3, 4))
                                        + model_mat_risk[..., list(self.infection_index[risk_ind].values())].sum((0, 2, 3, 4, 5))
                                        + model_mat_risk[..., self.natural_immunity_index[risk_ind]].sum((0, 2, 3, 4)))

        genotype_frequency_normal = np.where(total_for_normalization_normal != 0,
                                            (genotype_frequency_normal / total_for_normalization_normal) * 100,
                                            0)
        
        genotype_frequency_general = np.where(model_mat_risk.sum((0, 2, 3, 4)) != 0,
                                            (genotype_frequency_general / model_mat_risk.sum((0, 2, 3, 4))) * 100,
                                            0)
        return genotype_frequency_normal,genotype_frequency_CIN, genotype_frequency_general
        
    
    def calc_cancer_related_results(self, risk_ind):
        last_precancer_state = [self.precancer_index[risk_ind][i][-1] for i in range(self.num_subgroup_type)]    
        
        new_cancer_cases = self.model_change[risk_ind][..., last_precancer_state, self.cancer_index[risk_ind][0][0]].sum((2,3,4))
        new_cancer_diagnosis = self.model_change[risk_ind][..., self.cancer_index[risk_ind][0], self.cancer_index[risk_ind][1]].sum((2,3,4))
        
        # special cases when calculating cancer
        cancer_mortality_rate = config.get(CANCER_MORTALITY_RATE_KEY)
        natural_mortality_rate = config.get(NATURAL_MORTALITY_RATE_KEY)

        rate_proportion = np.array(self.rate_data[risk_ind][cancer_mortality_rate].values)
        rate_total = rate_proportion + np.array(self.rate_data[risk_ind][natural_mortality_rate].values)[:,None]
        rate_proportion/=rate_total
        new_cancer_deaths_mod = self.model_change[risk_ind][..., self.cancer_index[risk_ind][1], self.death_state[self.modeling_disease_id]].sum((2,3))
        new_cancer_deaths = (new_cancer_deaths_mod * rate_proportion).sum(-1)
        new_natural_deaths = (new_cancer_deaths_mod * (1-rate_proportion)).sum(-1)
        return new_cancer_cases, new_cancer_diagnosis, new_cancer_deaths, new_natural_deaths
    
    
    def calc_total_infected(self, risk_ind):
        return self.model_matrix[risk_ind][...,self.overall_infection_index[risk_ind]].sum((-3,-2,-1))
        
        
    
    
    