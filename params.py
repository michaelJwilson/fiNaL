import  yaml


As             = 2.e-9
ns             = 0.965
r              = 0.0

zmax           = 7.0
dz             = 0.05

H0             = 67.6           # [km/s/Mpc].                                                                                                                                                                                               
h              = H0 / 100.

c              = 2.9979 * 1.e5  # [km/s].                                                                                                                                                                                                   
omk            = 0.0
ombh2          = 0.022
omch2          = 0.122
mnu            = 0.000
tau            = 0.06

Om             = (omch2 + ombh2) / h /h


if __name__ == '__main__':
    vv             = [x for x in dir() if not x.startswith('_')]
    vv.remove('yaml')

    ll             = locals()

    defined        = {k: ll[k] for k in vv}
    
    with open('dat/params.yaml', 'w') as ff:
        yaml.dump(defined, ff, default_flow_style=False)
