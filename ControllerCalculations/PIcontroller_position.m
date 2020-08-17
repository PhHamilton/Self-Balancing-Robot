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
G = n*km / (s*(R*((J*n^2 + J2)*s + dm*n^2) + n^2*km*kemf));

figure()
step(G)

w1 = 8;
w2 = w1;

xi1 = 1.2;
xi2 = xi1;

Am = s^2 + 2*xi1*w1*s + w1^2;
A0 = s^2 + 2*xi2*w2*s + w2^2;

b = km*n/(J*n^2 + J2);
a = n^2*(dm*R + km*kemf)/(J*n^2 + J2);

r0 = (2*(xi1*w1 + xi2*w2)*R*(J*n^2 + J2)-n^2*(dm*R + km*kemf))/(R*(J*n^2 + J2));
s0 = (w1*w2)^2*(J*n^2 + J2)/(km*n);
s1 = 2*(xi1*w1*w2^2  + xi2*w2*w1^2)*R*(J*n^2 + J2)/(km*n);
s2 = ((w1^2 + w2^2 + 4*xi1*w1*xi2*w2)*(R*(J*n^2 + J2)) - n^2*(dm*R + km*kemf)*r0)/(km*n);


A = s^2 + n^2/(J*n^2+ J2)*(dm*R + km*kemf)*s;
B = km*n/(J*n^2+ J2);
S = s2*s^2 + s1*s + s0;
R = s*(s+r0);

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
ylabel('Position [rad]');
title('Step comparison');
legend('Output Feedback', 'Error Feedback', 'location', 'southeast');

subplot(2,1,2);
pzmap(minreal(Gc_OF))

print("Plots/Position_StepNPole.eps");

# Discrete:
Gff = c2d(T/R, dt, 'tustin');
Gc  = c2d(S/R, dt, 'tustin');
