import numpy as np

class Halo:
    '''
    Skeletron class for SSHalo and NFWHalo. The two share the same code except
    for the different density profiles and phase space position retrival
    methods.
    '''

    def __init__(self, s=1.5, beta=1., size=1e2, z_ssb=1., **kwargs):
        '''
        Parameters
        ----------
            s : float
                Accretion rate, s<6 for this model
            beta : float
                Dimensionless symmetron force strength
            size : float
                Turn-around radius today / (vacuum compton wavelength)
            z_ssb : float
                Redshift of SSB, z>0
        '''
        # Define concentration based on s, c=rs/R
        self.s = s

        self.f0 = np.pi*np.pi/8.
        self.f1 = 8./3. * (1.+z_ssb)**3. * beta**2. / size**2.
        self.a0 = size**2./2
        self.b0 = 3*np.pi*np.pi*np.pi / 4. / (1.+z_ssb)**3.

        self.init_halo(**kwargs)

    def _orbit_integrator(self, pos_0, ts, symmetron=None, full=False):
        from scipy.integrate import odeint

        #**STIFF**
        big_number = 1e20
        width = 1e-5

        if(symmetron is None):
            def f(X, t):
                # X[0] is lambda, X[1] is lambda'
                dx0 = X[1]
                dx1 = -self.f0 * t**(2.*self.s/3.)/(X[0]**2.) * self.M(X[0]/self.Lambda(t))
                dx1 += big_number / (1 + np.exp(X[0]/width))
                return [dx0, dx1]
        else:
            def f(X, t):
                # X[0] is lambda, X[1] is lambda'
                dx0 = X[1]
                dx1 = -self.f0 * t**(2.*self.s/3.)/(X[0]**2.) * self.M(X[0]/self.Lambda(t)) \
                    - self.f1 * symmetron.force(X[0]/self.Lambda(t), t)/self.Lambda(t)
                dx1 += big_number / (1 + np.exp(X[0]/width))
                return [dx0, dx1]

        X = odeint(f, pos_0, ts, hmax=0.01)

        if(full):
            return X

        return X[-1, 0], X[-1, 1]


    def Lambda(self, t):
        return t**(2./3. + 2.*self.s/9.)

    def a(self, t):
        return self.a0 * t**(4./3. + 4./9. * self.s)

    def b(self, t):
        return self.b0 * t**(-2.)

    def _pos(self, t_collapse, symmetron=None, verbose=False):
        # return phase space position for t_collapse, with or without symmetron.

        t_collapse = np.atleast_1d(t_collapse)
        n = t_collapse.size

        result = np.zeros((n, 2))

        for i in xrange(n):
            if(t_collapse[i]==1.):
                result[i, 0], result[i, 1] = [1., 0.]
                continue

            ts = np.linspace(t_collapse[i], 1, 2e5)
            pos_0 = [self.Lambda(t_collapse[i]), 0.]
            result[i, 0], result[i, 1] = self._orbit_integrator(pos_0, ts, symmetron)
            if(verbose):
                print i, "-th step done"

        return result

class NFWHalo(Halo):
    def init_halo(self):
        s = self.s
        if(s>6):
            raise ValueError("s>6 is not supported for NFW model.")

        const = s/(1.+s/3.)
        c_array = np.linspace(1e-3, 2, 1000)
        func = c_array/self.fNFW(c_array) * self.fNFW_p(c_array)

        self.c = c_array[np.abs(func - const).argmin()]

    def fNFW(self, xi):
        return np.log(1.+xi) - xi/(1.+xi)

    def fNFW_p(self, xi):
        return xi/(1.+xi)**2.

    def M(self, y):
        return self.fNFW(y*self.c)/self.fNFW(self.c)

    def P(self, y):
        return self.fNFW_p(y*self.c)*self.c / self.fNFW(self.c)/4./np.pi/y**2.

    def pos(self, t_collapse, symmetron=None, verbose=False):
        # return phase space position for t_collapse, with or without symmetron.
        # use _orbit_integrator from super.
        return self._pos(t_collapse, symmetron, verbose)

class SSHalo(Halo):
    def init_halo(self, load_dir="cache", extend='const'):
        from scipy.interpolate import interp1d

        y = np.load(load_dir+"/y.npy")
        P = np.load(load_dir+"/P.npy")
        M = np.load(load_dir+"/M.npy")
        ts = np.load(load_dir+"/ts.npy")
        l = np.load(load_dir+"/l.npy")
        self.s = np.load(load_dir+"/s.npy")

        self.M_interp = interp1d(y, M, fill_value='extrapolate')
        if(self.s==1.5) and (extend == 'const'):
            # Outside or range, assume P constant
            y_plus = np.linspace(1, 100, int(1e4))
            P_plus = P[-1]*np.ones(int(1e4))
            self.P_interp = interp1d(np.concatenate((y, y_plus)), np.concatenate((P, P_plus)), fill_value='extrapolate')
        else:
            self.P_interp = interp1d(y, P, fill_value='extrapolate')


        self._lambda = interp1d(ts, l[:, 0], fill_value='extrapolate')
        self._lambda_p = interp1d(ts, l[:, 1], fill_value='extrapolate')

    def M(self, y):
        return self.M_interp(y)

    def P(self, y):
        return self.P_interp(y)

    def pos(self, t_collapse, symmetron=None, verbose=False, load=True):
        if(load and (symmetron is None) ):
            t_collapse = np.atleast_1d(t_collapse)
            n = t_collapse.size

            result = np.zeros((n, 2))
            for i in xrange(n):
                result[i, 0] = self._lambda(1./t_collapse[i])/self.Lambda(1./t_collapse[i])
                result[i, 1] = self._lambda_p(1./t_collapse[i])/self.Lambda(1./t_collapse[i])/t_collapse[i]
            return result
        else:
            return self._pos(t_collapse, symmetron, verbose)

class symmetron:
    def __init__(self, halo, ts=np.linspace(0.1, 1, 10), y=np.linspace(1e-3, 2, 2e2), EPS=1e-7, omega=0.9):
        # ys must be linearly spaced
        from scipy.interpolate import interp1d
        from fsolver import symmetron_integrator

        ts = np.array(ts)
        y = np.array(y)
        self.ts = ts
        self.y = y
        self.chichip = [None]*len(ts)
        self.phi = [None]*len(ts)
        self.min = [None]*len(ts)
        self.res = [None]*len(ts)

        i=0
        for t in ts:
            print i
            P = halo.P(y)
            phi = np.ones(len(P))
            res = np.zeros(len(P))
            a = halo.a(t)
            b = halo.b(t)

            symmetron_integrator.solver(phi, res, P, y, a, b, EPS, omega, len(P))
            self.res[i] = interp1d(y, res, fill_value=0, bounds_error=False)
            self.phi[i] = interp1d(y, phi, fill_value=0, bounds_error=False)
            self.chichip[i] = interp1d(y, phi*np.gradient(phi, y[2]-y[1]), fill_value=0, bounds_error=False)
            min = np.sqrt(1.-b*P)
            min[~np.isfinite(min)] = 0
            self.min[i] = interp1d(y, min, fill_value=0, bounds_error=False)
            i+=1

    def force(self, y, t):
        ts = self.ts

        if(t < ts.min()):
            return 0

        indices = np.arange(len(ts))

        if(t>1):
            min_idx = -1
            max_idx = -1
        else:
            min_idx = indices[ts<=t][-1]
            max_idx = indices[ts>=t][0]
            delta_t = ts[max_idx] - ts[min_idx]

        if(min_idx != max_idx):
            return \
                ((t - ts[min_idx])*self.chichip[min_idx](y) + (ts[max_idx]-t)*self.chichip[max_idx](y))/delta_t
        else:
            return self.chichip[min_idx](y)
