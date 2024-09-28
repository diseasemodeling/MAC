import numpy as np
import pandas as pd
import pdb
import random
class PartnershipAssignment():
    def __init__(self, num_samples = 1000,
                       time_unit = 1,
                       modification = None):
        # define parameters
        self.num_risk = 3
        self.age_lb = np.asarray([13,18,25,30,35,40,45])
        self.age_ub = np.asarray([17,24,29,34,39,44,65])
        self.age_group_interval = self.age_ub - self.age_lb + 1
        self.num_age_group = len(self.age_lb)
        self.num_samples = num_samples
        self.max_power = 7
        self.num_ptnrs_max = [0]
        self.num_ptnrs_min = [0]
        self.time_unit = time_unit
        self.prop_concurrency =np.asarray([0,0,0])# np.asarray([0.09823667,0.1307616, 0.1394358])
        self.rate_concurrency =np.asarray([0.78,0.78,1.55])
        self.modification = modification 
        
        
    def run(self):
        self.calc_num_ptnrs_range()
        self.read_data()
        self.calc_num_ptnrs_change()
        self.assign_ptnrs()
        self.calc_num_ptnrs_actual()
        self.output_final_results()

    # read data 
    def read_data(self):
        # read CDF of lifetime partnership initiation
        self.num_ptnr_init = np.zeros((self.num_risk, self.num_age_group, self.num_degree))
        df = pd.ExcelFile("lifetime_ptnr_CDF.xlsx")
        for i in range(self.num_risk):
            self.num_ptnr_init[i] =  df.parse(sheet_name = "{0}".format(i), header=None).values
    
    # calculate number of partners 
    def calc_num_ptnrs_range(self):
        for i in range(self.max_power+1):
            self.num_ptnrs_max.append(2**i)
        self.num_ptnrs_max = np.asarray(self.num_ptnrs_max)
        self.num_degree = self.num_ptnrs_max.shape[0]
        for j in range(1,self.num_degree):
            self.num_ptnrs_min.append(self.num_ptnrs_max[j-1]+1)
        self.num_ptnrs_min = np.asarray(self.num_ptnrs_min)

    # calculate number of parnters change per unit of time
    def calc_num_ptnrs_change(self, write = True):
        # calculate number of lifetime partnership initiation within each age group
        self.num_ptnrs_init_ub = np.zeros_like(self.num_ptnr_init)
        self.num_ptnrs_init_lb = np.zeros_like(self.num_ptnr_init)
        self.num_ptnrs_init_actual = (self.num_ptnrs_max+self.num_ptnrs_min)/2
        if self.modification is not None: 
            print("with modification")
            self.num_ptnrs_init_actual = (self.num_ptnrs_max+self.num_ptnrs_min)/2/self.modification
            
        j = 0     
        ## using average number
        while j < self.num_age_group:
            if j == 0: 
                self.num_ptnrs_init_ub[:, j] = np.ceil(self.num_ptnr_init[:,j] * self.num_ptnrs_init_actual[np.newaxis,:])
            else:
                self.num_ptnrs_init_ub[:, j] = np.ceil(self.num_ptnr_init[:,j] * self.num_ptnrs_init_actual[np.newaxis,:]) - np.ceil(self.num_ptnr_init[:,j-1] * self.num_ptnrs_init_actual[np.newaxis,:])
                
            j += 1 
        
        self.avg_num_ptnrs_init = self.num_ptnrs_init_ub /self.age_group_interval[np.newaxis, :,np.newaxis]/ self.time_unit
        
        if write == True:
            with pd.ExcelWriter('contact_rate_change.xlsx') as writer:
                for i in range(self.num_risk):
                    pd.DataFrame(self.avg_num_ptnrs_init[i]).to_excel(writer, sheet_name='{0}'.format(i))

    # randomly assign the time step at which the lifetime partnership initiates
    def assign_ptnrs(self):
        self.num_ptnrs = np.zeros((self.num_samples, self.num_risk,self.time_unit*(self.age_ub[-1] - self.age_lb[0] +1),  self.num_degree))
        size =  self.num_ptnrs_init_ub # if you sum(0), you will find it equals to average total number of lifetime partners 

        for i in range(self.num_samples):
            for j in range(self.num_risk):
                for d in range(self.num_degree):
                    for a in range(self.num_age_group):
                        # number of ages by time unit 
                        num_time_steps = self.age_group_interval[a] * self.time_unit 
                        selected_ind = np.random.choice(num_time_steps, round(size[j,a,d]), replace = True)
                        start = (self.age_lb[a] - self.age_lb[0])* self.time_unit 
                        selected_ind += start
                        
                        for ind in selected_ind: 
                            self.num_ptnrs[i,j,ind,d] += 1
        
        
    # fill out all zero elements along age group after the first nonzero elements
    def calc_num_ptnrs_actual(self, concurrency = True):
        self.num_ptnrs_actual = np.copy(self.num_ptnrs)
        i = 0
        while i < self.num_samples:
            j = 0
            while j < self.num_risk:
                d = 0
                while d < self.num_degree:
                    a = 0
                    found = False
                    # find the first nonzero elements and a is the index 
                    while a < self.num_ptnrs.shape[2] and found == False:
                        if self.num_ptnrs[i,j,a,d] > 0:
                            found = True
                        a += 1
                    self.num_ptnrs_actual[i,j,a:,d][self.num_ptnrs_actual[i,j,a:,d] == 0] = 1 
                    d += 1
                j += 1
            i += 1
        if concurrency is True:
            print("with concurrency")
            self.update_concurrency()
        else:
            print("no concurrency")
    
    def output_final_results(self, write = True):
        avg_num_ptnrs = np.zeros((self.num_samples, self.num_risk, self.num_age_group, self.num_degree))

        for a in range(self.num_age_group):
            start, end = (self.age_lb[a] - self.age_lb[0])*self.time_unit, (self.age_ub[a] - self.age_lb[0] + 1)*self.time_unit
            avg_num_ptnrs[...,a,:] = np.mean(self.num_ptnrs_actual[...,start:end,:],axis = 2)
        self.avg_num_ptnrs = np.mean(avg_num_ptnrs, axis = 0)

        if write == True:
            with pd.ExcelWriter('avg_num_ptnrs.xlsx') as writer:
                for i in range(self.num_risk):
                    pd.DataFrame(self.avg_num_ptnrs[i]).to_excel(writer, sheet_name='{0}'.format(i))
                

    def update_concurrency(self):
        nonzeros_element = np.where(self.num_ptnrs > 0)
        total_nonzeros_element = len(nonzeros_element[0])
        
        self.duration_ptnrship = np.zeros_like(self.num_ptnrs)
        
        for length in range(total_nonzeros_element):
            i = nonzeros_element[0][length]
            j = nonzeros_element[1][length]
            a = nonzeros_element[2][length]
            d = nonzeros_element[3][length]
            
            size = int(self.num_ptnrs[i,j,a,d])
            is_concurrency = np.random.uniform(size = size) < self.prop_concurrency[j]
        
            if np.any(is_concurrency):
                len_concurrency = is_concurrency * np.round(np.random.exponential(1/self.rate_concurrency[j], size = size) * self.time_unit) 
                for m in len_concurrency:
                    if m != 0:
                        self.num_ptnrs_actual[i,j,a:(a+int(m)),d] += 1        
        
    # find age group the selected age belongs to 
    def find_age_index(self, age):
        i = 0
        found = False
        age_index = -1
        while i < len(self.age_ub) and found == False:
            if age <= self.age_ub[i]:
                age_index = i
                found = True
            i += 1
            
        ### just in case 
        if age > self.age_ub[-1]: 
            age_index = -1

        elif age < self.age_lb[0]: 
            age_index = 0

        return age_index

    def find_next_nonzero_elements(self, arr, ind):
        a = ind + 1
        while a < arr.shape[0]:
            if arr[a] > 0: return a
            a += 1
        return -1  


if __name__ == "__main__":
    pa = PartnershipAssignment(time_unit=2, num_samples=1000,modification =2.5)
    pa.run()
    

    