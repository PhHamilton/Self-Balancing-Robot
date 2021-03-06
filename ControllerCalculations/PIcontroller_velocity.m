clear all; close all;

updateFrequency = 20;
dt = 1/updateFrequency;

# Motor Parameters
R    = 1;
n    = 517;

noloadVoltage = 6;
noloadSpeed = 210;
stallTorque = 10e-2;
stallCurrent = 3.2;

km   = stallTorque/(stallCurrent);
kemf = 30*noloadVoltage/(pi*n*noloadSpeed);
J    = 2.5e-6;#4e-6;
J2   = 0;

dm   = 3e-6;

# System dynamics
s = tf('s');
G = n*km / (R*((J*n^2 + J2)*s + dm*n^2) + n^2*km*kemf);

# Controller:
  #Pole Placement
    w0 = 5;
    w1 = w0;

Am = s + w0;
A0 = s + w1;

# Calculation of parameters
s0 = w0*w1*R*(J*n^2 + J2)/(n*km);
s1 = ((w0 + w1)*R*(J*n^2 + J2) - (dm*R + km*kemf)*n^2)/(n*km);

#A = R*((J*n^2 + J2)*s + dm*n^2) + n^2*km*kemf;
A = s + n^2/(J*n^2+ J2)*(dm*R + km*kemf);
B = km*n/(J*n^2+ J2);

R = s;
S = s1*s + s0;


t0 = dcgain(Am)/B;
T  = t0*A0;

# Output Feedback
  Gc_OF = minreal(B*T/(A*R + B*S));

# Error Feedback:
  F = S/R;
  Gc_EF = minreal(F*G/(1 + F*G));

figure()
subplot(2,1,1);
step(Gc_OF, Gc_EF)
ylabel('Velocity [rad/s]');
title('Step comparison');
legend('Output Feedback', 'Error Feedback', 'location', 'southeast');

subplot(2,1,2);
pzmap(Gc_OF)

print("Plots/Velocity_StepNPole.eps");
# Control signal: u = T/R * r - S/R * y

### Discrete System ###
Gff = c2d(T/R, dt, 'tustin');
Gc  = c2d(S/R, dt, 'tustin');
