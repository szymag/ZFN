class Temperature(object):
    coefficients = {'c': (1., 0., 272.15), 'f': (1.8, -273.15, 32.), 'r': (1.8, 0., 0.)}

    def __init__(self, **kwargs):
        try:
            name, value = kwargs.popitem()
        except KeyError:
            name, value = 'k', 0.
        if kwargs or name not in 'kcfr':
            kwargs[name] = value
            raise TypeError('nieprawidłowe parametry %r' % kwargs)
        setattr(self, name, float(value))

    def __getattr__(self, name):
        try:
            eq = self.coefficients[name]
        except KeyError:
            raise AttributeError(name)
        return (self.k + eq[1]) * eq[0] + eq[2]

    def __setattr__(self, name, value):
        if name in self.coefficients:
            eq = self.coefficients[name]
            self.k = (value - eq[2]) / eq[0] - eq[1]
        elif name == 'k':
            object.__setattr__(self, name, value)
        else:
            raise AttributeError(name)

    def __str__(self):
        return '%s K' % self.k

    def __repr__(self):
        return 'Temperature(k=%r)' % self.k


t = Temperature(f=70)
print(t.f)
