# -*- coding: utf-8 -*-
"""
Izhikevich single neuron memory model 
Created on Tue Nov 14 20:37:26 2017

@author: Zhongxi

According to the article Simple Model of Spiking Neurons (Eugene M. Izhikevich)
There is certain combination of paramters that can make the single Izhikevich neuron able to memorize some states.
That is, with a predetermined fixed potentiation, we can "turn on" the neuron by appling a short stimulus pulse
and turn off the neuron with a negative pulse

This code manually finds such parameter combinations.

Attention: there is overshooting: 
        1) shutting down very high stimulus current can actually turn off the circuit
        2) shutting down very high inhibition can also turn on the neuron (!)
"""


from brian2 import *
import matplotlib.pyplot as plt
start_scope()

# %%    Parameters
N_neuron = 1
a = 0.1
b = 0.26
c = -65 # for future exploration: I'd rather keep this value fixed for biological reasons
d = 2   # for future exploration: I'd rather keep this value fixed for biological reasons

duration = 1500 # ms
# %%     Neuron model and monitors
neuron_eqs = '''
dv/dt=(0.04*v**2+5*v+140-u+I_input(t))/ms: 1
du/dt=(a*(b*v-u))/ms : 1
a : 1
b : 1
c : 1
d : 1
'''
Neuron_pop = NeuronGroup(N_neuron, neuron_eqs, method='euler',threshold = 'v>30', reset ='v=c \n u=u+d')
state_mon = StateMonitor(Neuron_pop,['v','u'],record=True)
spike_mon = SpikeMonitor(Neuron_pop)

Neuron_pop.a = a
Neuron_pop.b = b
Neuron_pop.c = c
Neuron_pop.d = d
# %%    Test
I_potentiation = 0.2    
t0 = 100                    # steady state time
I1 = 0.15                   # stimulus (in addition to potentiation)
t1 = 50                     # stimulus time
I2 = -I_potentiation        # inhibition
t2 = 50 # ms                # inhibition time 

vec0 = np.arange(t0)*0                  # wait to steady state
vec1 = I_potentiation + np.ones(150)*0   # potentiation
vec2 = I_potentiation + np.ones(t1)*I1  # stimulation
vec3 = I_potentiation + np.ones(200)*0  # observe
vec4 = I_potentiation + np.ones(t2)*I2  # inhibition
vec5 = I_potentiation + np.ones(duration)*0 # observe
I_array =  np.concatenate((vec0,vec1,vec2,vec3,vec4,vec5))
t_array =  np.arange(len(I_array))
I_input = TimedArray(I_array, dt=1*ms)
store()

run(duration*ms)
# %%
fig, ax1 = plt.subplots()
ax1.plot(state_mon.t/ms,state_mon.v[0],'b');
ax1.set_xlabel('time (ms)')
# Make the y-axis label, ticks and tick labels match the line color.
ax1.set_ylabel('v', color='b')
ax1.tick_params('y', colors='b')

ax2 = ax1.twinx()
ax2.plot(t_array*second,I_array,'r');
ax2.set_ylabel('I', color='r')
ax2.tick_params('y', colors='r')

fig.tight_layout()
plt.show()
xlim(100,1500)

figure()
plot(state_mon.t/ms,state_mon.u[0],'b');