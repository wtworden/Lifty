
from mpmath import iv
from math import ceil
from math import log


def convert(string_or_int):
    """convert a real interval from the string form '1.23244?' (used by Sage), 
        to an mp interval iv.mpf([1.23243,1.23245]). convert a complex interval from the string 
        form '1.232? + 2.34?*I' to the interval form iv.mpc([1.231,1.233],[2.33,2.35]).
    """

    c = string_or_int

    if isinstance(c,int) or isinstance(c,sage.rings.integer.Integer):
        return iv.mpc(int(c))
    else:
        assert isinstance(c,str)
        if c[0] != '-':
            c = '+' + c
        c = c.replace(' ','')
        c = c.replace('-','*-1*')
        c = c.replace('+','*1*')
        c = c.split('*')
        c = c[1:]

    
        if len(c) == 2:
            im = 0
            re = c[1]
            re_sign = int(c[0])
            if re.find('?') != -1:
                re = re.replace('?','')
                decs = len(re)-re.find('.')-1
                re_a = str((int(re.replace('.',''))-1))[:-decs]+'.'+str((int(re.replace('.',''))-1))[-decs:]
                re_b = str((int(re.replace('.',''))+1))[:-decs]+'.'+str((int(re.replace('.',''))+1))[-decs:]
                re = re_sign*iv.mpf([re_a,re_b])
            else:
                re = re_sign*iv.mpf(re)
            return iv.mpc(re,im)
    
        elif len(c) == 3:
            re = 0
            im = c[1]
            im_sign = int(c[0])
            if im.find('?') != -1:
                im = im.replace('?','')
                decs = len(im)-im.find('.')-1
                im_a = str((int(im.replace('.',''))-1))[:-decs]+'.'+str((int(im.replace('.',''))-1))[-decs:]
                im_b = str((int(im.replace('.',''))+1))[:-decs]+'.'+str((int(im.replace('.',''))+1))[-decs:]
                im = im_sign*iv.mpf([im_a,im_b])
            else:
                im = im_sign*iv.mpf(im)
            return iv.mpc(re,im)
    
        elif len(c) == 5:
            re = c[1]
            im = c[3]
            re_sign = int(c[0])
            im_sign = int(c[2])
            if re.find('?') != -1:
                re = re.replace('?','')
                decs = len(re)-re.find('.')-1
                re_a = str((int(re.replace('.',''))-1))[:-decs]+'.'+str((int(re.replace('.',''))-1))[-decs:]
                re_b = str((int(re.replace('.',''))+1))[:-decs]+'.'+str((int(re.replace('.',''))+1))[-decs:]
                re = re_sign*iv.mpf([re_a,re_b])
            else:
                re = re_sign*iv.mpf(re)
    
            if im.find('?') != -1:
                im = im.replace('?','')
                decs = len(im)-im.find('.')-1
                im_a = str((int(im.replace('.',''))-1))[:-decs]+'.'+str((int(im.replace('.',''))-1))[-decs:]
                im_b = str((int(im.replace('.',''))+1))[:-decs]+'.'+str((int(im.replace('.',''))+1))[-decs:]
                im = im_sign*iv.mpf([im_a,im_b])
            else:
                im = im_sign*iv.mpf(im)
    
            return iv.mpc(re,im)






class Polynomial:
    def __init__(self, coefficients):
        """ A Polynomial is specified by a list (or a tuple) of coefficients [a_n, a_{n-1}, ..., a_0]. 
            Each complex coefficient can be entered as follows:
            '1.23? + 4.56?*I'  (the question mark indicates that the last digits is correct up to +/- 1)
            '1.23 + 4.56*I'    (coefficient is known exactly)
            5                  (coefficient is an integer) Note: always enter decimals as strings
            '1 - 4.56?*I'      (real part known exactly, imaginary part in interval [4.55,4.57])
            '1.23' or '1.23?'  (real coefficient, exact or interval)
            '4.56*I' or '4.56?*I'  (real part is 0, exact or interval)

            For example, the co-rabbit polynomial can be entered as:

            >>> Polynomial([1,0,'-.12256? - .74486?*I'])
        """
        self.coefficients = []
        for c in coefficients:
            if not isinstance(c,iv.mpc):
                c = convert(c)
            self.coefficients.append(c)

        self.degree = len(self.coefficients)-1
        self._roots = None
    
    def coeff(self, k):
        return self.coefficients[self.degree-k]

    def evaluate(self,z):
        if not isinstance(z,iv.mpc):
            z = convert(z)
        result = iv.mpc(0)
        for i in range(self.degree+1):
            result += self.coeff(i)*(z**i)
        return result

    def iterate(self,z,k):
        """ evaluate P^k(z)
        """
        result = self.evaluate(z)
        for i in range(k-1):
            result = self.evaluate(result)
        return result

    def postcritical(self):
        pass

    def critial_pts(self):
        pass

    def derivative(self):
        Q_coeffs = []
        P_coeffs = self.coefficients
        for i in range(self.degree):
            n = self.degree - i
            Q_coeffs.append(n*P_coeffs[i])
        return Polynomial(Q_coeffs)

    def newton(self,z):
        return z - (self.evaluate(z))/(self.derivative().evaluate(z))


    def roots(self, epsilon=None):
        if self._roots == None:
            self._roots = self._find_roots(epsilon)
            return self._roots
        elif epsilon == None:
            return self._roots
        else:
            self._roots = self._find_roots(epsilon)
            return self._roots

    def _find_roots(self, epsilon):
        eps = epsilon
        d = self.degree

        if d == 1:
            return [-self.coeff(0)/self.coeff(1)]

        if d == 2:
            a, b, c = self.coefficients
            root1 = (-b + iv.sqrt(b**2-4*a*c))/(2*a)
            root2 = (-b - iv.sqrt(b**2-4*a*c))/(2*a)
            return [root1, root2]

        s = int(ceil(0.26632*log(d)))
        N = int(ceil(8.32547*d*log(d)))
        R_fl = 1 + sqrt(2)
        R = 1 + iv.sqrt(2)
        K = int(ceil(d*log(R_fl/eps)))
        Sd_next = []
        eps = iv.mpf(str(eps))

        for i in range(1,s+1):
            r = R*iv.mpf(((d-1)/d))**((2*i-1)/(4*s))
            for j in range(N):
                v = 2*iv.pi*j/N
                Sd_next.append(r*iv.exp(iv.mpc(0,v)))

        approx_roots = []

        while len(Sd_next)>0 and len(approx_roots)<d:  
            Sd_current = Sd_next
            Sd_next = []
            for z in Sd_current:
                z_old = z
                dispatched = False
                for i in range(K):
                    z_new = self.newton(z_old)
                    if iv.absmax(z_new-z_old) < eps/d:
                        root_re = iv.mpf([z_new.real.a-eps,z_new.real.b+eps])
                        root_im = iv.mpf([z_new.imag.a-eps,z_new.imag.b+eps])
                        new_root = iv.mpc(root_re,root_im)
                        for root in approx_roots:
                            if root.overlap(new_root):
                                dispatched = True
                                break
                        if not dispatched:
                            approx_roots.append(new_root)
                            dispatched = True
                        break
                    else:
                        if iv.isinf(iv.absmax(z_old)):
                            #raise IntervalError('Need to increase precision of polynomial coefficients, or decrease required root precision by increasing epsilon.')
                            dispatched=True
                            break
                        else:
                            z_old = z_new
                if not dispatched:
                    Sd_next.append(z_new)
        repeat_roots = []
        if len(approx_roots) < self.degree:
            Q = self.derivative()
            deriv_roots = Q.roots(epsilon)
            overlaps = dict([(deriv_root,[]) for deriv_root in deriv_roots])
            for deriv_root in deriv_roots:
                matched = False
                for root in approx_roots:
                    if deriv_root.overlap(root):
                        overlaps[deriv_root].append(root)
                        matched = True
                if not matched:
                    if self.evaluate(deriv_root).overlap(iv.mpc([0,0])):
                        repeat_roots.append(deriv_root)
            missing_repetitions = []
            for i in range(len(repeat_roots)):
                if len([repeat_roots[j] for j in range(i) if repeat_roots[i].overlap(repeat_roots[j])]) == 0:
                    missing_repetitions.append(repeat_roots[i])

            repeat_roots = repeat_roots + missing_repetitions

            for deriv_root in deriv_roots:
                if len(overlaps[deriv_root]) == 1:
                    repeat_roots.append(deriv_root)
                elif len(overlaps[deriv_root]) > 1:
                    raise IntervalError('Need to increase epsilon (which might require increases precision of coefficients, and iv.dps')
        all_roots = approx_roots + repeat_roots
        assert len(all_roots) == self.degree
        for root in all_roots:
            assert self.evaluate(root).overlap(iv.mpc([0,0]))

        return all_roots


    def postcrit_set(self):
        


class IntervalError(Exception):
    def __init__(self,msg):
        if msg:
            self.message = msg
        else:
            self.message = None

    def __str__(self):
        return 'IntervalError: {0}'.format(self.message)







