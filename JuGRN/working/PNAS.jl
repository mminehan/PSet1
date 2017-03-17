# Graph is just promoter concentration as a fcn of rxn time
# Transcription with Promoter lag
# R + P -k1/k-1-> RP_c -k2/k-2-> RP_o
# We know (assuming [R] >> [P]_total, SS for RP_c, and equil lies far to right)

using PyPlot

P_total = 10.0^-3 # arbitrary units

k_obs = 3.0*10.0^-3

n = 50000 # number of data points
t_beg = 0 # seconds, start time of simulation
t_end = 10*60 # seconds, end time of simulation

t = collect(t_beg:(t_end-t_beg)/n:t_end)

RP_o = P_total * (1-e.^-(k_obs*t))
RP_o_inc = P_total

max_transcription_rate = 60.0 #nt/s, from DataDictionary

product_conc = max_transcription_rate*RP_o.*t
product_conc_inc = max_transcription_rate*RP_o_inc.*t

plot(t,product_conc)
plot(t,product_conc_inc)

xlabel("Time (s)")
ylabel("[mRNA] ")


# 2C

L_t = 100 # bp, characteristic gene length
L_t1 = 10 # bp, subject gene length
L_t2 = 100
L_t3 = 1000

t_rate1 = max_transcription_rate*(L_t/L_t1)
t_rate2 = max_transcription_rate*(L_t/L_t2)
t_rate3 = max_transcription_rate*(L_t/L_t3)

figure()

mRNA_conc1 = t_rate1*RP_o.*t
mRNA_conc2 = t_rate2*RP_o.*t
mRNA_conc3 = t_rate3*RP_o.*t

plot(t,mRNA_conc1)
plot(t,mRNA_conc2)
plot(t,mRNA_conc3)

xlabel("Time (s)")
ylabel("[mRNA] ")
