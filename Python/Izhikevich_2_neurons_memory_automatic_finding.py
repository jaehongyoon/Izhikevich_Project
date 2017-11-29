# -*- coding: utf-8 -*-
"""
Izhikevich single neuron memory model 
Created on Tue Nov 14 20:37:26 2017

@author: Zhongxi

According to the article Simple Model of Spiking Neurons (Eugene M. Izhikevich)
There is certain combination of paramters that can make the single Izhikevich neuron able to memorize some states.
That is, with a predetermined fixed potentiation, we can "turn on" the neuron by appling a short stimulus pulse
and turn off the neuron with a negative pulse

This code finds parameter combinations of a and b.

Attention: there is overshooting: 
        1) shutting down very high stimulus current can actually turn off the circuit
        2) shutting down very high inhibition can also turn on the neuron (!)
"""


from brian2 import *
import matplotlib.pyplot as plt
start_scope()

#global spike_mon, Neuron_pop

def score(spike_recording):
    score_list = np.arange(len(spike_recording))
    for i in range(len(spike_recording)):
        spike_t = spike_recording[i]/ms
        n_spike_A = len([num for num in spike_t if num>100 and num<250])
        n_spike_B = len([num for num in spike_t if num>350 and num<500])
        n_spike_C = len([num for num in spike_t if num>600 and num<750])
        score_list[i] = (n_spike_B/(n_spike_A+1))* (n_spike_B/(n_spike_C+1))
    return score_list

def evolve(gene_A,gene_B,gene_scores,retain=0.2, cross=0.2,random_select=0.05):
    # population = list of parameters 
    sort_idx = np.argsort(gene_scores)[::-1];
    gene_scores = gene_scores[sort_idx];
    gene_A = gene_A[sort_idx]*1.
    gene_B = gene_B[sort_idx]*1.
    N_population = len(gene_scores)
    #   1) retain good parents
    N_retain = int(retain*N_population)
    new_gene_A = gene_A[:N_retain-1]
    new_gene_B = gene_B[:N_retain-1]
    #   2) crossover to generate children
    N_children = int(N_population*cross)
    while len(new_gene_A) < N_children + N_retain:
        female = np.random.randint(0, N_population-1)
        male   = np.random.randint(0, N_population-1)
        if male != female:
            new_gene_A=np.append(new_gene_A,gene_A[male])
            new_gene_B=np.append(new_gene_B,gene_B[female])
    #   3) random selection
    N_randsel = int(N_population*random_select)
    while len(new_gene_A) < N_children + N_retain+N_randsel:
        rand_idx = np.random.randint(N_retain,N_population-1)
        new_gene_A=np.append(new_gene_A,gene_A[rand_idx])
        new_gene_B=np.append(new_gene_B,gene_B[rand_idx])
    #   4) mutation (Gaussian)
    std_A = np.std(new_gene_A[:N_retain-1])
    std_B = np.std(new_gene_B[:N_retain-1])
    new_gene_A=np.append(new_gene_A, clip(np.random.randn(N_population-len(new_gene_A))*std_A+new_gene_A[0], 0.01,1))
    new_gene_B=np.append(new_gene_B, clip(np.random.randn(N_population-len(new_gene_B))*std_B+new_gene_B[0], 0.01,1))
    return new_gene_A,new_gene_B,gene_scores,gene_A,gene_B

def show_evolution(gene_A,gene_B,gene_scores):
    for i in range(len(gene_scores)):
        if gene_scores[i]>1.5:
            plot(gene_A[i],gene_B[i],'.',color=np.array([1,0,0])*(1-(gene_scores[i])/40.))
    return
# %%    Parameters
N_neuron = 1
a = 0.1
b = 0.2
c = -65 # for future exploration: I'd rather keep this value fixed for biological reasons
d = 2  # for future exploration: I'd rather keep this value fixed for biological reasons
a_lim = [.1, 1.]
b_lim = [.1, 1.]
duration = 3000 # ms
Esyn = 0
tau_syn_d = 0.5*ms
tau_syn_r = 0.5*ms
w = 3000    #   For resonator: 1000~1500 is a decent range. others are possible. below <500 is impossible
            #   For fast spiking setting, below 2500 is not oscillating, above 2500 is self-oscillating. No controllable memory here.
            
# %% To-do's
'''
1) score function
2) general genetic algorithm
...
aim: analog memory 
'''
# %%     Neuron model and monitors
neuron_eqs_1 = '''
dv/dt=(0.04*v**2+5*v+140-u+I_input(t)+Isyn)/ms: 1
du/dt=(a*(b*v-u))/ms : 1

Isyn = g * (Esyn-v) : 1
dg/dt = ( -g/(tau_syn_d) +z*Hz ) : 1
dz/dt = ( -z/(tau_syn_r) ) : 1

a : 1
b : 1
c : 1
d : 1
'''
neuron_eqs_2 = '''
dv/dt=(0.04*v**2+5*v+140-u+Isyn)/ms: 1
du/dt=(a*(b*v-u))/ms : 1

Isyn = g * (Esyn-v) : 1
dg/dt = ( -g/(tau_syn_d) + z*Hz) : 1
dz/dt = ( -z/(tau_syn_r) ) : 1

a : 1
b : 1
c : 1
d : 1
'''
Neuron_pop1 = NeuronGroup(N_neuron, neuron_eqs_1, method='euler',threshold = 'v>30', reset ='v=c \n u=u+d') 
Neuron_pop2 = NeuronGroup(N_neuron, neuron_eqs_2, method='euler',threshold = 'v>30', reset ='v=c \n u=u+d') 
Syn_12 = Synapses(Neuron_pop1,Neuron_pop2,method='euler',on_pre='z_post += w')
Syn_21 = Synapses(Neuron_pop2,Neuron_pop1,method='euler',on_pre='z_post += w')
Syn_12.connect('i==j')
Syn_21.connect('i==j')

state_mon1 = StateMonitor(Neuron_pop1,['v','u'],record=True)
state_mon2 = StateMonitor(Neuron_pop2,['v','u','g'],record=True)
spike_mon = SpikeMonitor(Neuron_pop2)
Neuron_pop1.a = a
Neuron_pop1.b = b
Neuron_pop1.c = c
Neuron_pop1.d = d
Neuron_pop2.a = a
Neuron_pop2.b = b
Neuron_pop2.c = c
Neuron_pop2.d = d
# %%    Test stimulation
I_potentiation = 0.2    
t0 = 100                    # steady state time
t_stim = 100                     # stimulus time
I_depolarize = -I_potentiation        # inhibition
I_stim = np.array([0.15, 0.3, 0.45])*30   # stimulation list
t_depolarize = 50

vec0 = np.arange(t0)*0                  # wait to steady state
vec1 = I_potentiation + np.ones(150)*0   # potentiation
I_array = np.concatenate((vec0,vec1))
for idx in range(len(I_stim)):
    I_array = np.concatenate((I_array,  I_potentiation + np.ones(t_stim)*I_stim[idx]  ))
    I_array = np.concatenate((I_array,  I_potentiation + np.ones(200)*0  ))
    I_array = np.concatenate((I_array,  I_potentiation + np.ones(t_depolarize)*I_depolarize  ))
    I_array = np.concatenate((I_array,  I_potentiation + np.ones(200)*0  ))

I_input = TimedArray(I_array, dt=1*ms)
store()

run(duration*ms)

figure(1)
plot(state_mon1.t,state_mon1.v[0])
figure(2)
plot(state_mon2.t,state_mon2.v[0])

## %% Evolution
#iteration = 0
#Neuron_pop.a[0] = 0.1
#Neuron_pop.b[0] = 0.26
#Neuron_pop.a[1:] = np.random.rand(N_neuron-1)+0.01
#Neuron_pop.b[1:] = np.random.rand(N_neuron-1)+0.01
#next_gene_A = Neuron_pop.a*1.
#next_gene_B = Neuron_pop.b*1.
#scores = np.zeros(N_neuron)
#while np.min(scores)==0:
#    restore()
#    Neuron_pop.a = next_gene_A
#    Neuron_pop.b = next_gene_B
#    run(duration*ms)
#    scores = score(spike_mon.spike_trains())
#    next_gene_A,next_gene_B,sorted_scores,sorted_gene_A,sorted_gene_B = \
#    evolve(Neuron_pop.a,Neuron_pop.b,scores,retain=0.2, cross=0.2,random_select=0.05)
##    print('iter=%d a=%.2f b=%.2f'%(iteration,para_A[0],para_B[0]))
#    show_evolution(sorted_gene_A,sorted_gene_B,sorted_scores)
#    pause(1*ms)
#    iteration += 1
## %% debugging
##fig, ax1 = plt.subplots()
##ax1.plot(state_mon.t/ms,state_mon.v[0],'b');
##ax1.set_xlabel('time (ms)')
### Make the y-axis label, ticks and tick labels match the line color.
##ax1.set_ylabel('v', color='b')
##ax1.tick_params('y', colors='b')
##
##ax2 = ax1.twinx()
##ax2.plot(t_array*second,I_array,'r');
##ax2.set_ylabel('I', color='r')
##ax2.tick_params('y', colors='r')
##
##fig.tight_layout()
##plt.show()
##xlim(100,1500)
