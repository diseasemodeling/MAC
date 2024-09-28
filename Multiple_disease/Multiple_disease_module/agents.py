import numpy as np
from .model import Model
import pickle
import pdb
from numpy.random import default_rng



def find_age_index(age, age_lb, age_ub):
    age_index = 0
    for i in range(len(age_lb)):
        if age > age_lb[i] and age <= age_ub[i]:
            age_index = i
            break
    if age > age_ub[-1]:
        age_index = len(age_ub) - 1
    # if age < age_lb[0]:
    #     age_index = 0
    # else: age_index = np.searchsorted(age_ub, age, side='right')
    return age_index

def find_degree_bin(degree, degree_dist_bin_all):
    return np.searchsorted(degree_dist_bin_all, degree)

def adjust_rounding_matrix(input, target, rng):
    # Initial rounding
    input_round = np.round(input)
    
    # Get a mask of non-zero values
    mask_nonzero = input != 0
    
    # Calculate the residual difference between target and the sum of the rounded values along the last axis
    res = target - input_round.sum(-1)
    
    # Adjust values: 
    # Generate random offsets for each entry, sort them, and select as many as needed to adjust the total
    # This step can produce both negative and positive adjustments
    input_offset = np.array((-rng.random(input_round.shape) * \
                             mask_nonzero).argsort(-1) < np.abs(res)[..., None],\
                             dtype=np.float64) * np.sign(res)[..., None]
    
    # Apply adjustments
    input_round += input_offset
    
    # Ensure all values are positive integers
    # If any value goes negative, adjust it and ensure the matrix still sums up to the target
    while np.any(input_round < 0):
        for i in np.argwhere(input_round < 0):
            input_round[tuple(i)] = 0
            available_to_borrow = np.where(input_round[tuple(i[:-1])] > 1)
            if available_to_borrow[0].size > 0:
                idx_to_reduce = available_to_borrow[0][0]
                coord_to_reduce = tuple(list(i[:-1]) + [idx_to_reduce])
                input_round[coord_to_reduce] -= 1
                input_round[tuple(i)] += 1
    
    return input_round

class Agents(Model):
    def __init__(self, time_unit = 2, 
                        mode = "netlogo",
                        increased_risk = 0,
                        initial_year = 2006,
                        dryrun_step = 70,
                        simulation_step = 10,
                        model_screening = False,
                        model_vaccine = False,
                        plot = False,
                        random_seed = 1076):
        super().__init__(time_unit = time_unit, 
                        mode = mode,
                        increased_risk = increased_risk,
                        initial_year = initial_year,
                        dryrun_step = dryrun_step,
                        simulation_step = simulation_step,
                        model_screening = model_screening,
                        model_vaccine = model_vaccine, 
                        plot= plot,
                        random_seed= random_seed)
        
        # define variables from disease
        self.base_disease_id = self.disease.base_disease_id
        self.modeling_disease_id = self.disease.modeling_disease_id
        self.num_disease_state = self.disease.num_disease_state
        self.death_state = self.disease.death_state
        self.susceptible_index = self.disease.susceptible_index
        self.rng = default_rng(self.random_seed)
        # initialize variables
        self.initialize_variables()
        
    def initialize_variables(self):
        self.init_agent_matrix = {n:np.zeros((self.num_disease_state[self.base_disease_id][n], 
                                            self.num_age_group, 
                                            self.num_degree_bin, 
                                            self.num_jurisdictions))
                                  for n in range(self.num_risk_group)}
        self.state_assignment = {n:np.zeros((self.num_age_group, 
                                           self.num_degree_bin, 
                                           self.num_jurisdictions, 
                                           self.num_disease_state[self.modeling_disease_id][n]))
                                 for n in range(self.num_risk_group)} 
    def collect_agents_info(self, 
                            age, 
                            degree, 
                            sexInd,
                            jur_ind, 
                            base_disease_state):
        age_ind = find_age_index(age, age_lb = self.age_lb, age_ub = self.age_ub)
        degree_ind = find_degree_bin(degree, degree_dist_bin_all = self.degree_dist_bin_all)
        self.init_agent_matrix[sexInd][base_disease_state, age_ind, degree_ind, jur_ind] += 1
    
    def determine_multi_disease_transition(self, 
                                           risk_ind, 
                                           age_ind, 
                                           degree_ind, 
                                           modeling_state, 
                                           base_state,
                                           jur_ind):
        if modeling_state == -1:
            total_agents = self.state_assignment[risk_ind][age_ind, degree_ind, jur_ind]
        else:
            total_agents = self.model_change[risk_ind][base_state, age_ind, degree_ind, jur_ind, modeling_state]

        non_zero_indices = np.where(total_agents>0)
        if len(non_zero_indices[0]) > 0:
            if len(non_zero_indices[0]) > 1: 
                self.rng.shuffle(non_zero_indices[0])
            inflow_values = total_agents[non_zero_indices].astype(int)
            # if np.any(total_agents[non_zero_indices] < inflow_values): print("total",total_agents[non_zero_indices], "round", inflow_values)
            # if np.any(np.asarray(list(inflow_values))<0): print(list(non_zero_indices[0]), list(inflow_values), base_state, age_ind, degree_ind, jur_ind, modeling_state)
            return list(non_zero_indices[0]), list(inflow_values)

        else: 
            return [],[]

    
    def disease_state_agents_assignment(self, risk_ind):
        self.model_matrix[risk_ind][self.nonagents_state, ...,self.death_state[self.modeling_disease_id][risk_ind]] = 0 # manually reset number of death to 0
        # manually set place where number of people < 1 to 0
        np.maximum(self.model_matrix[risk_ind], 0, out=self.model_matrix[risk_ind])
        disease_dist = np.nan_to_num(self.model_matrix[risk_ind][self.nonagents_state]/self.model_matrix[risk_ind][self.nonagents_state].sum(-1)[...,None])
        disease_dist = np.nan_to_num(disease_dist/disease_dist.sum(-1)[...,None])
        self.state_assignment[risk_ind] = self.init_agent_matrix[risk_ind][...,None].sum(0) * disease_dist 
        self.state_assignment[risk_ind] = adjust_rounding_matrix(self.state_assignment[risk_ind], self.init_agent_matrix[risk_ind].sum(0), self.rng)
        
    
    def disease_state_nonagents_assignment(self,risk_ind, num_nodes = None):
        model = np.delete(self.model_matrix[risk_ind][self.nonagents_state], self.death_state[self.modeling_disease_id][risk_ind], axis = -1)
        model_sum = model.sum(-1)[...,None]
        disease_dist = np.nan_to_num(model/model_sum)
        index = np.delete(np.arange(self.num_disease_state[self.modeling_disease_id][risk_ind]), self.death_state[self.modeling_disease_id][risk_ind])
        if num_nodes is None: num_nodes = self.num_nodes
        num_nodes =  np.asarray(num_nodes)
        num_nodes_expanded = num_nodes[risk_ind, None, None, None, None]  
        degree_dist_expanded = self.degree_dist_by_age_degree[risk_ind,...,None, None] 
        self.model_matrix[risk_ind][self.nonagents_state][...,index] = disease_dist * num_nodes_expanded * degree_dist_expanded
            
    def disease_state_assignment(self, num_nodes = None):
        for risk_ind in range(self.num_risk_group):
            self.disease_state_nonagents_assignment(risk_ind, num_nodes)
            self.disease_state_agents_assignment(risk_ind)

    # update multi-disease state in the compartment for agents  
    def update_multi_disease_state(self, 
                                   base_disease_state, 
                                   age, 
                                   degree, 
                                   sexInd, 
                                   jur_ind, 
                                   modeling_state, 
                                   disease_id, 
                                   new_state):
        # Find indices for age and degree
        age_index = find_age_index(age, age_lb=self.age_lb, age_ub=self.age_ub)
        degree_index = find_degree_bin(degree, degree_dist_bin_all=self.degree_dist_bin_all)
        
        # Access the relevant matrix cell
        temp_mat = self.model_matrix[sexInd][base_disease_state, age_index, degree_index, jur_ind, modeling_state]

        # Check if the matrix value needs updating and is not in a default state
        if temp_mat < 1 and modeling_state != -1:pass

        # Decrement the matrix value if not in a default state
        if modeling_state != -1 and temp_mat >= 1:
            self.model_matrix[sexInd][base_disease_state, age_index, degree_index, jur_ind, modeling_state] -= 1         
            # if disease_id == self.modeling_disease_id: self.netlogo_change[sexInd][base_disease_state, age_index, degree_index, jur_ind, modeling_state, modeling_state] -= 1
        
        # Update the model matrix based on disease id
        if disease_id == self.base_disease_id:
            self.model_matrix[sexInd][new_state, age_index, degree_index, jur_ind, modeling_state] += 1
        elif disease_id == self.modeling_disease_id:
            # Special handling for modeling disease state
            if modeling_state == -1 or temp_mat >= 1:
                self.model_matrix[sexInd][base_disease_state, age_index, degree_index, jur_ind, new_state] += 1
                # self.netlogo_change[sexInd][base_disease_state, age_index, degree_index, jur_ind, modeling_state, new_state] += 1
        else:
            raise ValueError("Invalid disease id.")

        # Verify population after updating disease state
        self.check_population("update_disease_state")

       
    ### update aging or removing aging out         
    def update_aging(self, 
                     sexInd, 
                     base_disease_state, 
                     modeling_state, 
                     degree, 
                     jur_ind, 
                     old_age_ind, 
                     new_age_ind = None):
        degree_ind = find_degree_bin(degree, self.degree_dist_bin_all)
        old_age_value = self.model_matrix[sexInd][base_disease_state, old_age_ind, degree_ind, jur_ind, modeling_state]
        if old_age_value >= 1: 
            self.model_matrix[sexInd][base_disease_state, old_age_ind, degree_ind, jur_ind, modeling_state] = max(0, old_age_value - 1)
            
            if new_age_ind is not None:  
                self.model_matrix[sexInd][base_disease_state, new_age_ind, degree_ind, jur_ind, modeling_state] += 1
            else:
                last_age_index = -1
                self.model_matrix[sexInd][self.nonagents_state, last_age_index, degree_ind, jur_ind, modeling_state] += 1
        else:pass
            # print("update_aging",old_age_value, base_disease_state, old_age_ind, degree_ind, jur_ind, modeling_state, new_age_ind)
        self.check_population("update_aging")
        
    # determine disease state of selected neigbors and update compartment matrix  
    def determine_newly_infected_disease_state(self, sexInd, age, degree, jur_ind, death):
        age_group = find_age_index(age,self.age_lb, self.age_ub)
        degree_bin = find_degree_bin(degree, self.degree_dist_bin_all)
       
        model_mat_sex = self.model_matrix[sexInd]
        temp_mat = np.delete(model_mat_sex[self.nonagents_state], self.death_state[self.modeling_disease_id][sexInd], axis=-1)  # Select from alive & compartment people
        temp_mat_slice = temp_mat[age_group, degree_bin, jur_ind]
        
        # Default state selection
        selected_state = self.susceptible_index[sexInd][0]

        # Calculate total population
        total_pop = temp_mat_slice.sum()
        
        # Handle state selection and modification
        if total_pop < 1 or np.all(temp_mat_slice < 1):
            # print("total population < 1",temp_mat_slice, total_pop, sexInd, age, degree)
            temp_mat_slice[temp_mat_slice < 0] = 0
        else:
            if np.any(temp_mat_slice < 0): 
                # print(temp_mat) 
                temp_mat_slice[temp_mat_slice < 0] = 0
            # Normalize and randomly select state
            state_list = np.where(temp_mat_slice > 1)
            p_val = temp_mat_slice[state_list]/temp_mat_slice[state_list].sum()
            selected_state = self.rng.choice(state_list[0], p=p_val)
            model_mat_sex[self.nonagents_state, age_group, degree_bin, jur_ind, selected_state] -= 1

        # Increase count for infected state
        if death:
            newly_infected_ind = -1
            model_mat_sex[newly_infected_ind, age_group, degree_bin, jur_ind, selected_state] += 1
        else:
            newly_infected_ind = 1
            model_mat_sex[newly_infected_ind, age_group, degree_bin, jur_ind, selected_state] += 1
        
        return selected_state
    

    def write_HIV_results(self, list1):
        self.output.new_HIV_inf.append(list1[0])
        self.output.new_HIV_diag.append(list1[1])
