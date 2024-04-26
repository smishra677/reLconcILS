import gym
from gym import spaces
import numpy as np
import copy
import sys
sys.path.append("./reconcILS/")
sys.path.append('./')
import utils.Tree as Tree
import copy
import utils.Tally as Tally
import argparse
import utils.ILS as ILS
import utils.readWrite as readWrite
import pickle
import gc
import time
from reconcILS import *
from concurrent.futures import ThreadPoolExecutor, TimeoutError
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


tollerance= 0.001  
discount= 1.0 -1E-3 
alpha=0.01          

class ReLconcILS(gym.Env):
    metadata = {'render.modes': ['console']}

    def __init__(self, sp_string,gene_tree):
        super(ReLconcILS, self).__init__()
        self.id=0
        self.gene_string=gene_tree
        self.sp_string=sp_string
        self.red= readWrite.readWrite()
        self.species_tree = self.red.parse(sp_string)
        self.gene_tree =self.red.parse(gene_tree)
        self.list_mode=self.get_node()
        self.action_space = spaces.MultiDiscrete(3)
        self.observation_space = spaces.Box(low=0, high=255, shape=(1,), dtype=np.uint8)
        self.state_dic={}
 


    def get_node(self):
        count=[]
        stack=[self.species_tree]
        while stack:
            current_node = stack.pop()
            if current_node:
                count+=[current_node]
                if current_node.rightChild:
                    stack.append(current_node.rightChild)
                if current_node.leftChild:
                    stack.append(current_node.leftChild)
        return count

   

    def step(self, action):
        if len(self.list_mode)==0:
            self.state_dic[self.id]=self.species_tree
            self.current_state = self.id
            self.id=self.id+1
            return self.current_state, -10, True, {}


        print(action)
  
        index= action[1]

        node=self.list_mode[index]
        print('node',node.to_newick(),self.species_tree.to_newick())


        operation_type = action[0]
        

        if operation_type == 0:
            self._apply_duplication(node)
        elif operation_type == 1:
            self._apply_loss(node)
        elif operation_type == 2:
            self.perform_nni(node)
        
        self.list_mode=self.get_node()
        

        if len(self.list_mode)==0:
            self.state_dic[self.id]=self.species_tree
            self.current_state = self.id
            self.id=self.id+1
            return self.current_state, -100, True, {}

        if self.species_tree:

            self.species_tree= self.red.parse(self.species_tree.to_newick())
        done = self._check_done()


        if not done:

            reward = -(10*self._count_discordant_branches(self.gene_tree,self.species_tree)+1)

        else:
            
            reward=-1
            try:
                reward += 100
            except:
                reward=-100

        self.state_dic[self.id]=self.species_tree
        self.current_state = self.id
        self.id=self.id+1
            
        
        info = {}
        
        return self.current_state, reward, done, info

    def _apply_duplication(self, sp):
        recon= reconcils()
        recon_left,recon_right = sp.deepcopy_double()


        sp.leftChild=recon_left
        sp.rightChild=recon_right

        recon.clearid(recon_left,'Left')
                                
                                
                                

        recon.clearid(recon_right,'right')
                                       
                                        
        sp.children+=[sp.leftChild ]
        sp.children+=[sp.rightChild]
        sp.taxa=''
        self.species_tree=self.red.parse(sp.to_newick())


        

                                    

    def _apply_loss(self, sp):
        if sp.parent and sp.isLeaf==None:
            print(sp.parent.to_newick())
            if sp.parent.leftChild==sp:
                sp.parent.leftChild=None
                self.species_tree=sp
            else:
                sp.parent.rightChild=None
                self.species_tree=sp
        else:
            del sp
            self.species_tree=None

        
        if self.species_tree != None and self.species_tree.to_newick()==';':
            self.species_tree=None


        
        



        

        


        

    def perform_nni(self, sp):
        if sp and sp.isLeaf==None and sp.leftChild and sp.rightChild:
            print(sp.taxa)
            lis_l,lis_r=None,None
            print(self.species_tree.to_newick())
            print(sp.to_newick())
            sp.label_internal()
            if sp.leftChild.isLeaf==None :
                lis_l =sp.NNI(self.species_tree,'Left')
            
            if sp.rightChild.isLeaf==None:
                lis_r =sp.NNI(self.species_tree,'Right')

            value_dic= {}
            value_list=[]
            if lis_l:
                value_list =[lis_l[0][1],lis_l[1][1]]
            
                value_dic[0] = self._count_discordant_branches(self.gene_tree,lis_l[0][1])
                value_dic[1] = self._count_discordant_branches(self.gene_tree,lis_l[1][1])
            
            if lis_r:
                if len(value_list)==0:
                    value_dic[0]= self._count_discordant_branches(self.gene_tree,lis_r[0][1])
                    value_dic[1]= self._count_discordant_branches(self.gene_tree,lis_r[1][1])
                else:
                    value_dic[3]= self._count_discordant_branches(self.gene_tree,lis_r[0][1])
                    value_dic[4]= self._count_discordant_branches(self.gene_tree,lis_r[1][1])


                value_list+=[lis_r[0][1],lis_r[1][1]]



            if len(value_list)!=0: 
                min_key = min(value_dic, key=value_dic.get)
                self.species_tree= value_list[min_key]


        

    

    def count_(self,sp):
        count=0
        stack=[sp]
        while stack:
            current_node = stack.pop()
            if current_node:

                count+=len(current_node.refTo)

                if current_node.rightChild:
                    stack.append(current_node.rightChild)
                if current_node.leftChild:
                    stack.append(current_node.leftChild)
        return count

        

    def _count_discordant_branches(self,sp,tr):
        sp.reset()
        tr.reset()
        tr.order_gene(sp)
        tr.label_internal()
        sp.label_internal()
        tr.map_gene(sp)
        return self.count_(sp)
        





    def _check_done(self):
        if self.species_tree:
            print(self.species_tree.to_newick(),self.gene_tree.to_newick())
            self.gene_tree.reset()
            self.species_tree.reset()
            self.species_tree.order_gene(self.gene_tree)
            self.species_tree.label_internal()
            self.gene_tree.label_internal()
            self.species_tree.map_gene(self.gene_tree)
            count=[]
            stack=[self.gene_tree]
            while stack:
                current_node = stack.pop()
                if current_node:
                    count.append(len(current_node.refTo))

                    if current_node.rightChild:
                        stack.append(current_node.rightChild)
                    if current_node.leftChild:
                        stack.append(current_node.leftChild)
            
            if len(count)==sum(count):
                return True
            elif self.species_tree==None:
                return True
            else:
                return False
        return True

    def reset(self):
        red= readWrite.readWrite()
        self.species_tree = red.parse(self.sp_string)
        self.gene_tree =red.parse(self.gene_string)
        self.current_state = self.species_tree
        self.id=0
        self.list_mode=self.get_node()
        return self.id,None



def epsilon_greedy(env,Q, state, epsilon, nA):
    if np.random.rand() < epsilon:
        rand1 =np.random.choice(nA)
        rand2 = np.random.choice(len(env.list_mode))

        while rand1*rand2>=len(env.list_mode):
            rand1 =np.random.choice(nA)
            rand2 = np.random.choice(len(env.list_mode))




        return [rand1,rand2]
    else: 
        rand= np.random.choice(np.flatnonzero(Q[state] == np.max(Q[state])))
        action= rand%3 
        branch = int(rand/3)
        while branch>=len(env.list_mode):
            print(branch,len(env.list_mode))
            rand= np.random.choice(np.flatnonzero(Q[state] == np.max(Q[state])))
            action= rand%3 
            branch = int(rand/3)
        return [action,branch]
    
visit_count={}



def update_Q(Q, state, action, reward,terminated, next_state, next_action, alpha, gamma):
    action=action[0]*action[1]
    next_action= next_action[0]*next_action[1]
    if terminated:
        Q[state, action] =Q[state, action] + alpha* reward

    else:
        Q[state, action] += alpha * (reward + gamma * Q[next_state, next_action] - Q[state, action])
      

def evaluate_policy(env, policy,gamma):
    total_rewards = []
    for i in range(5):
        state, prob = env.reset()
        print(i)

        episode_reward = 0
        isfinal = False
        epi=0
        while not isfinal and epi<5:
            
            rand= np.argmax(policy[state])
            action= rand%3 
            branch = int(rand/3)
            while branch>=len(env.list_mode):
                rand= np.argmax(policy[state])
                action= rand%3 
                branch = int(rand/3)
                print(rand,branch)
            action= [action,branch]
           
            state, reward, isfinal, _ = env.step(action)
            epi=epi+1
            episode_reward +=((gamma)**epi)*reward
        total_rewards.append(episode_reward)
        print('-------------------------------------------------')
    return sum(total_rewards) / 5

def update_val_Q(env,nA,Q):
    Q[:,nA*len(env.list_mode) :][Q[:,nA*len(env.list_mode) :] != -1e6] = -1e6

    return Q

def sarsa(env, nS, nA,rewards,ep, episodes, alpha=0.01, gamma=0.999):

    Q = np.zeros((nS, nA*20))
    Q = update_val_Q(env,nA,Q)
    data=[]

    for I in range(episodes):
        epsilon  = ep*(episodes-I-1)/(episodes-1)
        state, _ = env.reset()


        action = epsilon_greedy(env,Q, state, epsilon, nA)
        print(I,action,state)
        done = False
        step=0
        while not done:

           
            next_state, reward, terminated, _ = env.step(action)
            done=terminated
            if done:
                break
            print(I,action,next_state,step)

            
            next_action = epsilon_greedy(env,Q, next_state, epsilon, nA)

            update_Q(Q, state, action, reward,done, next_state, next_action, alpha, gamma)
            Q = update_val_Q(env,nA,Q)
            
            state, action = next_state, next_action

            
            #data.append
            step=step+1
            print(I)


        if I%2==0:
            print('evaluate')

            if I in rewards:
                rewards[I]+=[evaluate_policy(env,Q,gamma)]
            else:
                rewards[I]=[evaluate_policy(env,Q,gamma)]
        


        

    return Q, rewards

if __name__ == "__main__":
    env = ReLconcILS('(A,(B,C));','(A,A)')

    
    reward ={}
    for i in range(2):
        env = ReLconcILS('(A,(B,C));','(C,(A,B))')
        Q,reward =sarsa(env, 100, 3,reward,ep=0.5,episodes=100, alpha=0.01, gamma=0.999)
   
    env = ReLconcILS('(A,(B,C));','(B,(A,C))')
    evaluate_policy(env, Q,0.999)
    df = pd.DataFrame(reward)

    ax= df.boxplot()
    ax.set_title('sarsa')
    fig = ax.get_figure()
    fig.savefig('sarsa_reconcILS.png', dpi=300)
    plt.clf()

