import numpy as np
import matplotlib.pyplot as plt
from scipy.special import iv
from scipy.integrate import quad_vec


a = 14.054
b = 0.473
c_ = np.array([
    [-3/4, 23/64, -7/64],
    [-1/2, 17/16, -9/16]
])


def I_0(x):
    return iv(0, x)


def I_1(x):
    return iv(1, x)


def Q_G(r, kappa):
    return (3/4/np.pi/kappa)**(3/2)*np.exp(-3*r**2/4/kappa)


def Q_D(r, kappa):
    return Q_G(r, kappa)*(1 - 5*kappa/4 + 2*r**2 - 33*r**4/80/kappa)


def J_SY(kappa):
    return 112.04*kappa**2*np.exp(0.246/kappa - a*kappa)


def J_SYD(kappa):
    if kappa > 1/8:
        return J_SY(kappa)
    else:
        return Q_D(0, kappa)

def Q_I(r, kappa):
    c = 1 - (1 + (0.38*kappa**(-0.95))**(-5))**(-1/5)
    if kappa < 1/8:
        d = 1
    else:
        d = 1 - 1/(0.177/(kappa - 0.111) + 6.40*(kappa - 0.111)**0.783)
    sum = 0
    for i in (-1, 0):
        for j in (1, 2, 3):
            sum += c_[i + 1, j - 1]*kappa**i*r**(2*j)
    h = (1 - c*r**2)/(1 - r**2)
    t = r/(1 - b**2*r**2)
    arg = -d*kappa*a*(1 + b)*t
    return J_SYD(kappa) * \
        (h)**(5/2) * \
        np.exp(sum/(1 - r**2)) * \
        np.exp(arg*b*r) * \
        I_0(arg)


def d_Q_I_d_r(r, kappa):
    c = 1 - (1 + (0.38*kappa**(-0.95))**(-5))**(-1/5)
    if kappa < 1/8:
        d = 1
    else:
        d = 1 - 1/(0.177/(kappa - 0.111) + 6.40*(kappa - 0.111)**0.783)
    sum = 0
    d_sum_dr = 0
    for i in (-1, 0):
        for j in (1, 2, 3):
            sum += c_[i + 1, j - 1]*kappa**i*r**(2*j)
            d_sum_dr += 2*j*c_[i + 1, j - 1]*kappa**i*r**(2*j - 1)
    h = (1 - c*r**2)/(1 - r**2)
    t = r/(1 - b**2*r**2)
    arg = -d*kappa*a*(1 + b)*t
    return J_SYD(kappa) * \
        (h)**(5/2) * \
        np.exp(sum/(1 - r**2)) * \
        np.exp(arg*b*r) * (
            (
                5/2/h*(2*r*h/(1 - r**2) - 2*r*c/(1 - r**2)) +
                d_sum_dr/(1 - r**2) + 2*r*sum/(1 - r**2)**2 +
                2*b*arg *(
                    1 + r*b**2*t
                )
            ) * I_0(arg)
            + arg*(
                1/r + 2*b**2*t
            ) * I_1(arg)
        )


gamma = np.linspace(0, 1, 10000)[1:-1]
kappas = (1/7, 5/27, 5/19, 5/11, 5/3)
for kappa in kappas:
    P_eq = Q_I(gamma, kappa)
    plt.plot(gamma, P_eq, label=kappa)
plt.legend()
plt.show()
for kappa in kappas:
    P_eq = Q_I(gamma, kappa)
    g_eq = 4*np.pi*gamma**2*P_eq
    plt.plot(gamma, g_eq, label=kappa)
plt.legend()
plt.show()
for kappa in kappas:
    P_eq = Q_I(gamma, kappa)
    beta_Delta_psi = -np.log(P_eq/P_eq[0])
    plt.semilogy(gamma, beta_Delta_psi, label=kappa)
plt.legend()
plt.show()
for kappa in kappas:
    eta = -d_Q_I_d_r(gamma, kappa)/Q_I(gamma, kappa)
    plt.semilogy(gamma, eta, label=kappa)
    P_eq = Q_I(gamma, kappa)
    beta_Delta_psi = -np.log(P_eq/P_eq[0])
    eta_numerical = np.gradient(beta_Delta_psi)/np.gradient(gamma)
    plt.semilogy(gamma, eta_numerical, 'k:')
    Z, _ = quad_vec(
        lambda gamma:
            Q_I(gamma, kappa)*np.sinh(eta*gamma)*gamma/eta,
        1e-6, 1 - 1e-6
    )
    gamma_isotensional = np.gradient(np.log(Z))/np.gradient(eta)
    plt.semilogy(gamma_isotensional, eta, 'k--')
plt.legend()
plt.show()