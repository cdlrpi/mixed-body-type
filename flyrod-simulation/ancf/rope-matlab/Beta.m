function beta = Beta(e1,e2,e3,e4,e5,e6,e7,e8,l,rho,A,Y,II)
                    
g = 9.81;

beta1 = -(Y / l ^ 3 * (288 * A * e1 ^ 3 - 288 * A * e5 ^ 3 - 3360 * II * e5 + 3360 * II * e1 - 14 * A * l ^ 3 * e7 - 864 * A * e1 ^ 2 * e5 + 864 * A * e1 * e5 ^ 2 + 288 * A * e1 * e2 ^ 2 + 288 * A * e1 * e6 ^ 2 - 288 * A * e5 * e2 ^ 2 - 288 * A * e5 * e6 ^ 2 - A * l ^ 3 * e7 ^ 3 - A * l ^ 3 * e3 ^ 3 - 168 * A * e1 * l ^ 2 + 168 * A * e5 * l ^ 2 - 14 * A * l ^ 3 * e3 + 1680 * II * l * e7 + 1680 * II * l * e3 + 108 * A * e1 ^ 2 * l * e7 + 36 * A * e1 * l ^ 2 * e7 ^ 2 + 12 * A * e1 * l ^ 2 * e4 ^ 2 + 36 * A * e1 * l ^ 2 * e3 ^ 2 + 108 * A * e1 ^ 2 * l * e3 + 12 * A * e1 * l ^ 2 * e8 ^ 2 + 108 * A * e5 ^ 2 * l * e7 - 36 * A * e5 * l ^ 2 * e7 ^ 2 - 12 * A * e5 * l ^ 2 * e4 ^ 2 - 36 * A * e5 * l ^ 2 * e3 ^ 2 + 108 * A * e5 ^ 2 * l * e3 - 12 * A * e5 * l ^ 2 * e8 ^ 2 + A * l ^ 3 * e7 * e4 ^ 2 + 3 * A * l ^ 3 * e7 * e3 ^ 2 - A * l ^ 3 * e7 * e8 ^ 2 + 3 * A * l ^ 3 * e7 ^ 2 * e3 + 36 * A * l * e7 * e2 ^ 2 + 36 * A * l * e7 * e6 ^ 2 - A * l ^ 3 * e3 * e4 ^ 2 + A * l ^ 3 * e3 * e8 ^ 2 + 36 * A * l * e3 * e2 ^ 2 + 36 * A * l * e3 * e6 ^ 2 - 576 * A * e1 * e2 * e6 + 576 * A * e5 * e2 * e6 - 216 * A * e1 * l * e5 * e7 + 72 * A * e1 * l * e2 * e4 - 216 * A * e1 * l * e3 * e5 + 72 * A * e1 * l * e2 * e8 - 72 * A * e1 * l * e6 * e8 - 72 * A * e1 * l * e4 * e6 - 72 * A * e5 * l * e2 * e4 - 72 * A * e5 * l * e2 * e8 + 72 * A * e5 * l * e6 * e8 + 72 * A * e5 * l * e4 * e6 + 24 * A * l ^ 2 * e7 * e2 * e8 - 72 * A * l * e7 * e2 * e6 + 2 * A * l ^ 3 * e7 * e4 * e8 - 24 * A * l ^ 2 * e7 * e6 * e8 + 24 * A * l ^ 2 * e3 * e2 * e4 - 72 * A * l * e3 * e2 * e6 + 2 * A * l ^ 3 * e3 * e4 * e8 - 24 * A * l ^ 2 * e3 * e4 * e6)) / 0.280e3;
beta2 = -l * rho * A * g / 0.2e1 - Y / l ^ 3 * (-0.14e2 * A * l ^ 3 * e8 + 0.24e2 * A * l ^ 2 * e8 * e1 * e7 - 0.24e2 * A * l ^ 2 * e8 * e5 * e7 - 0.72e2 * A * l * e8 * e1 * e5 - 0.216e3 * A * l * e8 * e2 * e6 + 0.2e1 * A * l ^ 3 * e8 * e3 * e7 - 0.72e2 * A * l * e4 * e1 * e5 + 0.24e2 * A * l ^ 2 * e4 * e1 * e3 - 0.24e2 * A * l ^ 2 * e4 * e3 * e5 - 0.216e3 * A * l * e4 * e2 * e6 + 0.2e1 * A * l ^ 3 * e4 * e3 * e7 + 0.72e2 * A * e2 * l * e1 * e7 - 0.72e2 * A * e2 * l * e5 * e7 + 0.72e2 * A * e2 * l * e1 * e3 - 0.72e2 * A * e2 * l * e3 * e5 - 0.72e2 * A * e6 * l * e1 * e7 + 0.72e2 * A * e6 * l * e5 * e7 - 0.72e2 * A * e6 * l * e1 * e3 + 0.72e2 * A * e6 * l * e3 * e5 - A * l ^ 3 * e8 ^ 3 - A * l ^ 3 * e4 ^ 3 - 0.168e3 * A * e2 * l ^ 2 + 0.168e3 * A * e6 * l ^ 2 - 0.864e3 * A * e2 ^ 2 * e6 + 0.288e3 * A * e2 * e1 ^ 2 + 0.288e3 * A * e2 * e5 ^ 2 + 0.864e3 * A * e2 * e6 ^ 2 - 0.288e3 * A * e6 * e1 ^ 2 - 0.288e3 * A * e6 * e5 ^ 2 - 0.14e2 * A * l ^ 3 * e4 + 0.1680e4 * II * l * e4 + 0.1680e4 * II * l * e8 - 0.576e3 * A * e2 * e1 * e5 + 0.576e3 * A * e6 * e1 * e5 + 0.3e1 * A * l ^ 3 * e8 * e4 ^ 2 + A * l ^ 3 * e8 * e3 ^ 2 + 0.36e2 * A * l ^ 2 * e8 ^ 2 * e2 + 0.3e1 * A * l ^ 3 * e8 ^ 2 * e4 - 0.36e2 * A * l ^ 2 * e8 ^ 2 * e6 + 0.36e2 * A * l * e8 * e1 ^ 2 + 0.36e2 * A * l * e8 * e5 ^ 2 + 0.108e3 * A * l * e8 * e2 ^ 2 + 0.108e3 * A * l * e8 * e6 ^ 2 + A * l ^ 3 * e4 * e7 ^ 2 + 0.36e2 * A * l ^ 2 * e4 ^ 2 * e2 - A * l ^ 3 * e4 * e3 ^ 2 - 0.36e2 * A * l ^ 2 * e4 ^ 2 * e6 + 0.36e2 * A * l * e4 * e1 ^ 2 + 0.36e2 * A * l * e4 * e5 ^ 2 + 0.108e3 * A * l * e4 * e2 ^ 2 + 0.108e3 * A * l * e4 * e6 ^ 2 + 0.12e2 * A * e2 * l ^ 2 * e7 ^ 2 + 0.12e2 * A * e2 * l ^ 2 * e3 ^ 2 - 0.12e2 * A * e6 * l ^ 2 * e7 ^ 2 - 0.12e2 * A * e6 * l ^ 2 * e3 ^ 2 - A * l ^ 3 * e8 * e7 ^ 2 + 0.288e3 * A * e2 ^ 3 - 0.288e3 * A * e6 ^ 3 - 0.3360e4 * II * e6 + 0.3360e4 * II * e2) / 0.280e3;
beta3 = -(Y / l ^ 2 * (108 * A * e1 ^ 3 - 108 * A * e5 ^ 3 - 5040 * II * e5 + 5040 * II * e1 + 14 * A * l ^ 3 * e7 - 324 * A * e1 ^ 2 * e5 + 324 * A * e1 * e5 ^ 2 + 108 * A * e1 * e2 ^ 2 + 108 * A * e1 * e6 ^ 2 - 108 * A * e5 * e2 ^ 2 - 108 * A * e5 * e6 ^ 2 - 3 * A * l ^ 3 * e7 ^ 3 + 24 * A * l ^ 3 * e3 ^ 3 - 42 * A * e1 * l ^ 2 + 42 * A * e5 * l ^ 2 - 56 * A * l ^ 3 * e3 + 1680 * II * l * e7 + 3360 * II * l * e3 + 9 * A * e1 * l ^ 2 * e7 ^ 2 - 3 * A * e1 * l ^ 2 * e4 ^ 2 - 9 * A * e1 * l ^ 2 * e3 ^ 2 + 108 * A * e1 ^ 2 * l * e3 + 3 * A * e1 * l ^ 2 * e8 ^ 2 - 9 * A * e5 * l ^ 2 * e7 ^ 2 + 3 * A * e5 * l ^ 2 * e4 ^ 2 + 9 * A * e5 * l ^ 2 * e3 ^ 2 + 108 * A * e5 ^ 2 * l * e3 - 3 * A * e5 * l ^ 2 * e8 ^ 2 - 3 * A * l ^ 3 * e7 * e4 ^ 2 - 9 * A * l ^ 3 * e7 * e3 ^ 2 - 3 * A * l ^ 3 * e7 * e8 ^ 2 + 6 * A * l ^ 3 * e7 ^ 2 * e3 + 24 * A * l ^ 3 * e3 * e4 ^ 2 + 2 * A * l ^ 3 * e3 * e8 ^ 2 + 36 * A * l * e3 * e2 ^ 2 + 36 * A * l * e3 * e6 ^ 2 - 216 * A * e1 * e2 * e6 + 216 * A * e5 * e2 * e6 + 6 * A * e1 * l ^ 2 * e4 * e8 + 18 * A * e1 * l ^ 2 * e3 * e7 - 6 * A * e5 * l ^ 2 * e4 * e8 - 18 * A * e5 * l ^ 2 * e3 * e7 + 6 * A * l ^ 2 * e7 * e2 * e4 - 6 * A * l ^ 2 * e7 * e4 * e6 + 6 * A * l ^ 2 * e3 * e2 * e8 - 6 * A * l ^ 2 * e3 * e6 * e8 + 72 * A * e1 * l * e2 * e4 - 216 * A * e1 * l * e3 * e5 - 72 * A * e1 * l * e4 * e6 - 72 * A * e5 * l * e2 * e4 + 72 * A * e5 * l * e4 * e6 + 6 * A * l ^ 2 * e7 * e2 * e8 + 4 * A * l ^ 3 * e7 * e4 * e8 - 6 * A * l ^ 2 * e7 * e6 * e8 - 6 * A * l ^ 2 * e3 * e2 * e4 - 72 * A * l * e3 * e2 * e6 - 6 * A * l ^ 3 * e3 * e4 * e8 + 6 * A * l ^ 2 * e3 * e4 * e6)) / 0.840e3;
beta4 = -l ^ 2 * rho * A * g / 0.12e2 - Y / l ^ 2 * (0.14e2 * A * l ^ 3 * e8 + 0.6e1 * A * l ^ 2 * e8 * e1 * e7 - 0.6e1 * A * l ^ 2 * e8 * e5 * e7 + 0.4e1 * A * l ^ 3 * e8 * e3 * e7 - 0.72e2 * A * l * e4 * e1 * e5 - 0.6e1 * A * l ^ 2 * e4 * e1 * e3 + 0.6e1 * A * l ^ 2 * e4 * e3 * e5 - 0.216e3 * A * l * e4 * e2 * e6 - 0.6e1 * A * l ^ 3 * e4 * e3 * e7 + 0.72e2 * A * e2 * l * e1 * e3 - 0.72e2 * A * e2 * l * e3 * e5 - 0.72e2 * A * e6 * l * e1 * e3 + 0.72e2 * A * e6 * l * e3 * e5 + 0.18e2 * A * l ^ 2 * e8 * e2 * e4 + 0.6e1 * A * l ^ 2 * e8 * e1 * e3 - 0.6e1 * A * l ^ 2 * e8 * e3 * e5 - 0.18e2 * A * l ^ 2 * e8 * e4 * e6 + 0.6e1 * A * l ^ 2 * e4 * e1 * e7 - 0.6e1 * A * l ^ 2 * e4 * e5 * e7 + 0.6e1 * A * e2 * l ^ 2 * e3 * e7 - 0.6e1 * A * e6 * l ^ 2 * e3 * e7 - 0.3e1 * A * l ^ 3 * e8 ^ 3 + 0.24e2 * A * l ^ 3 * e4 ^ 3 - 0.42e2 * A * e2 * l ^ 2 + 0.42e2 * A * e6 * l ^ 2 - 0.324e3 * A * e2 ^ 2 * e6 + 0.108e3 * A * e2 * e1 ^ 2 + 0.108e3 * A * e2 * e5 ^ 2 + 0.324e3 * A * e2 * e6 ^ 2 - 0.108e3 * A * e6 * e1 ^ 2 - 0.108e3 * A * e6 * e5 ^ 2 - 0.56e2 * A * l ^ 3 * e4 + 0.3360e4 * II * l * e4 + 0.1680e4 * II * l * e8 - 0.216e3 * A * e2 * e1 * e5 + 0.216e3 * A * e6 * e1 * e5 - 0.9e1 * A * l ^ 3 * e8 * e4 ^ 2 - 0.3e1 * A * l ^ 3 * e8 * e3 ^ 2 + 0.9e1 * A * l ^ 2 * e8 ^ 2 * e2 + 0.6e1 * A * l ^ 3 * e8 ^ 2 * e4 - 0.9e1 * A * l ^ 2 * e8 ^ 2 * e6 + 0.2e1 * A * l ^ 3 * e4 * e7 ^ 2 - 0.9e1 * A * l ^ 2 * e4 ^ 2 * e2 + 0.24e2 * A * l ^ 3 * e4 * e3 ^ 2 + 0.9e1 * A * l ^ 2 * e4 ^ 2 * e6 + 0.36e2 * A * l * e4 * e1 ^ 2 + 0.36e2 * A * l * e4 * e5 ^ 2 + 0.108e3 * A * l * e4 * e2 ^ 2 + 0.108e3 * A * l * e4 * e6 ^ 2 + 0.3e1 * A * e2 * l ^ 2 * e7 ^ 2 - 0.3e1 * A * e2 * l ^ 2 * e3 ^ 2 - 0.3e1 * A * e6 * l ^ 2 * e7 ^ 2 + 0.3e1 * A * e6 * l ^ 2 * e3 ^ 2 - 0.3e1 * A * l ^ 3 * e8 * e7 ^ 2 + 0.108e3 * A * e2 ^ 3 - 0.108e3 * A * e6 ^ 3 - 0.5040e4 * II * e6 + 0.5040e4 * II * e2) / 0.840e3;
beta5 = (Y / l ^ 3 * (288 * A * e1 ^ 3 - 288 * A * e5 ^ 3 - 3360 * II * e5 + 3360 * II * e1 - 14 * A * l ^ 3 * e7 - 864 * A * e1 ^ 2 * e5 + 864 * A * e1 * e5 ^ 2 + 288 * A * e1 * e2 ^ 2 + 288 * A * e1 * e6 ^ 2 - 288 * A * e5 * e2 ^ 2 - 288 * A * e5 * e6 ^ 2 - A * l ^ 3 * e7 ^ 3 - A * l ^ 3 * e3 ^ 3 - 168 * A * e1 * l ^ 2 + 168 * A * e5 * l ^ 2 - 14 * A * l ^ 3 * e3 + 1680 * II * l * e7 + 1680 * II * l * e3 + 108 * A * e1 ^ 2 * l * e7 + 36 * A * e1 * l ^ 2 * e7 ^ 2 + 12 * A * e1 * l ^ 2 * e4 ^ 2 + 36 * A * e1 * l ^ 2 * e3 ^ 2 + 108 * A * e1 ^ 2 * l * e3 + 12 * A * e1 * l ^ 2 * e8 ^ 2 + 108 * A * e5 ^ 2 * l * e7 - 36 * A * e5 * l ^ 2 * e7 ^ 2 - 12 * A * e5 * l ^ 2 * e4 ^ 2 - 36 * A * e5 * l ^ 2 * e3 ^ 2 + 108 * A * e5 ^ 2 * l * e3 - 12 * A * e5 * l ^ 2 * e8 ^ 2 + A * l ^ 3 * e7 * e4 ^ 2 + 3 * A * l ^ 3 * e7 * e3 ^ 2 - A * l ^ 3 * e7 * e8 ^ 2 + 3 * A * l ^ 3 * e7 ^ 2 * e3 + 36 * A * l * e7 * e2 ^ 2 + 36 * A * l * e7 * e6 ^ 2 - A * l ^ 3 * e3 * e4 ^ 2 + A * l ^ 3 * e3 * e8 ^ 2 + 36 * A * l * e3 * e2 ^ 2 + 36 * A * l * e3 * e6 ^ 2 - 576 * A * e1 * e2 * e6 + 576 * A * e5 * e2 * e6 - 216 * A * e1 * l * e5 * e7 + 72 * A * e1 * l * e2 * e4 - 216 * A * e1 * l * e3 * e5 + 72 * A * e1 * l * e2 * e8 - 72 * A * e1 * l * e6 * e8 - 72 * A * e1 * l * e4 * e6 - 72 * A * e5 * l * e2 * e4 - 72 * A * e5 * l * e2 * e8 + 72 * A * e5 * l * e6 * e8 + 72 * A * e5 * l * e4 * e6 + 24 * A * l ^ 2 * e7 * e2 * e8 - 72 * A * l * e7 * e2 * e6 + 2 * A * l ^ 3 * e7 * e4 * e8 - 24 * A * l ^ 2 * e7 * e6 * e8 + 24 * A * l ^ 2 * e3 * e2 * e4 - 72 * A * l * e3 * e2 * e6 + 2 * A * l ^ 3 * e3 * e4 * e8 - 24 * A * l ^ 2 * e3 * e4 * e6)) / 0.280e3;
beta6 = -l * rho * A * g / 0.2e1 + Y / l ^ 3 * (-0.14e2 * A * l ^ 3 * e8 + 0.24e2 * A * l ^ 2 * e8 * e1 * e7 - 0.24e2 * A * l ^ 2 * e8 * e5 * e7 - 0.72e2 * A * l * e8 * e1 * e5 - 0.216e3 * A * l * e8 * e2 * e6 + 0.2e1 * A * l ^ 3 * e8 * e3 * e7 - 0.72e2 * A * l * e4 * e1 * e5 + 0.24e2 * A * l ^ 2 * e4 * e1 * e3 - 0.24e2 * A * l ^ 2 * e4 * e3 * e5 - 0.216e3 * A * l * e4 * e2 * e6 + 0.2e1 * A * l ^ 3 * e4 * e3 * e7 + 0.72e2 * A * e2 * l * e1 * e7 - 0.72e2 * A * e2 * l * e5 * e7 + 0.72e2 * A * e2 * l * e1 * e3 - 0.72e2 * A * e2 * l * e3 * e5 - 0.72e2 * A * e6 * l * e1 * e7 + 0.72e2 * A * e6 * l * e5 * e7 - 0.72e2 * A * e6 * l * e1 * e3 + 0.72e2 * A * e6 * l * e3 * e5 - A * l ^ 3 * e8 ^ 3 - A * l ^ 3 * e4 ^ 3 - 0.168e3 * A * e2 * l ^ 2 + 0.168e3 * A * e6 * l ^ 2 - 0.864e3 * A * e2 ^ 2 * e6 + 0.288e3 * A * e2 * e1 ^ 2 + 0.288e3 * A * e2 * e5 ^ 2 + 0.864e3 * A * e2 * e6 ^ 2 - 0.288e3 * A * e6 * e1 ^ 2 - 0.288e3 * A * e6 * e5 ^ 2 - 0.14e2 * A * l ^ 3 * e4 + 0.1680e4 * II * l * e4 + 0.1680e4 * II * l * e8 - 0.576e3 * A * e2 * e1 * e5 + 0.576e3 * A * e6 * e1 * e5 + 0.3e1 * A * l ^ 3 * e8 * e4 ^ 2 + A * l ^ 3 * e8 * e3 ^ 2 + 0.36e2 * A * l ^ 2 * e8 ^ 2 * e2 + 0.3e1 * A * l ^ 3 * e8 ^ 2 * e4 - 0.36e2 * A * l ^ 2 * e8 ^ 2 * e6 + 0.36e2 * A * l * e8 * e1 ^ 2 + 0.36e2 * A * l * e8 * e5 ^ 2 + 0.108e3 * A * l * e8 * e2 ^ 2 + 0.108e3 * A * l * e8 * e6 ^ 2 + A * l ^ 3 * e4 * e7 ^ 2 + 0.36e2 * A * l ^ 2 * e4 ^ 2 * e2 - A * l ^ 3 * e4 * e3 ^ 2 - 0.36e2 * A * l ^ 2 * e4 ^ 2 * e6 + 0.36e2 * A * l * e4 * e1 ^ 2 + 0.36e2 * A * l * e4 * e5 ^ 2 + 0.108e3 * A * l * e4 * e2 ^ 2 + 0.108e3 * A * l * e4 * e6 ^ 2 + 0.12e2 * A * e2 * l ^ 2 * e7 ^ 2 + 0.12e2 * A * e2 * l ^ 2 * e3 ^ 2 - 0.12e2 * A * e6 * l ^ 2 * e7 ^ 2 - 0.12e2 * A * e6 * l ^ 2 * e3 ^ 2 - A * l ^ 3 * e8 * e7 ^ 2 + 0.288e3 * A * e2 ^ 3 - 0.288e3 * A * e6 ^ 3 - 0.3360e4 * II * e6 + 0.3360e4 * II * e2) / 0.280e3;
beta7 = -(Y / l ^ 2 * (108 * A * e1 ^ 3 - 108 * A * e5 ^ 3 - 5040 * II * e5 + 5040 * II * e1 - 56 * A * l ^ 3 * e7 - 324 * A * e1 ^ 2 * e5 + 324 * A * e1 * e5 ^ 2 + 108 * A * e1 * e2 ^ 2 + 108 * A * e1 * e6 ^ 2 - 108 * A * e5 * e2 ^ 2 - 108 * A * e5 * e6 ^ 2 + 24 * A * l ^ 3 * e7 ^ 3 - 3 * A * l ^ 3 * e3 ^ 3 - 42 * A * e1 * l ^ 2 + 42 * A * e5 * l ^ 2 + 14 * A * l ^ 3 * e3 + 3360 * II * l * e7 + 1680 * II * l * e3 + 108 * A * e1 ^ 2 * l * e7 - 9 * A * e1 * l ^ 2 * e7 ^ 2 + 3 * A * e1 * l ^ 2 * e4 ^ 2 + 9 * A * e1 * l ^ 2 * e3 ^ 2 - 3 * A * e1 * l ^ 2 * e8 ^ 2 + 108 * A * e5 ^ 2 * l * e7 + 9 * A * e5 * l ^ 2 * e7 ^ 2 - 3 * A * e5 * l ^ 2 * e4 ^ 2 - 9 * A * e5 * l ^ 2 * e3 ^ 2 + 3 * A * e5 * l ^ 2 * e8 ^ 2 + 2 * A * l ^ 3 * e7 * e4 ^ 2 + 6 * A * l ^ 3 * e7 * e3 ^ 2 + 24 * A * l ^ 3 * e7 * e8 ^ 2 - 9 * A * l ^ 3 * e7 ^ 2 * e3 + 36 * A * l * e7 * e2 ^ 2 + 36 * A * l * e7 * e6 ^ 2 - 3 * A * l ^ 3 * e3 * e4 ^ 2 - 3 * A * l ^ 3 * e3 * e8 ^ 2 - 216 * A * e1 * e2 * e6 + 216 * A * e5 * e2 * e6 + 6 * A * e1 * l ^ 2 * e4 * e8 + 18 * A * e1 * l ^ 2 * e3 * e7 - 6 * A * e5 * l ^ 2 * e4 * e8 - 18 * A * e5 * l ^ 2 * e3 * e7 + 6 * A * l ^ 2 * e7 * e2 * e4 - 6 * A * l ^ 2 * e7 * e4 * e6 + 6 * A * l ^ 2 * e3 * e2 * e8 - 6 * A * l ^ 2 * e3 * e6 * e8 - 216 * A * e1 * l * e5 * e7 + 72 * A * e1 * l * e2 * e8 - 72 * A * e1 * l * e6 * e8 - 72 * A * e5 * l * e2 * e8 + 72 * A * e5 * l * e6 * e8 - 6 * A * l ^ 2 * e7 * e2 * e8 - 72 * A * l * e7 * e2 * e6 - 6 * A * l ^ 3 * e7 * e4 * e8 + 6 * A * l ^ 2 * e7 * e6 * e8 + 6 * A * l ^ 2 * e3 * e2 * e4 + 4 * A * l ^ 3 * e3 * e4 * e8 - 6 * A * l ^ 2 * e3 * e4 * e6)) / 0.840e3;
beta8 = l ^ 2 * rho * A * g / 0.12e2 - Y / l ^ 2 * (-0.56e2 * A * l ^ 3 * e8 - 0.6e1 * A * l ^ 2 * e8 * e1 * e7 + 0.6e1 * A * l ^ 2 * e8 * e5 * e7 - 0.72e2 * A * l * e8 * e1 * e5 - 0.216e3 * A * l * e8 * e2 * e6 - 0.6e1 * A * l ^ 3 * e8 * e3 * e7 + 0.6e1 * A * l ^ 2 * e4 * e1 * e3 - 0.6e1 * A * l ^ 2 * e4 * e3 * e5 + 0.4e1 * A * l ^ 3 * e4 * e3 * e7 + 0.72e2 * A * e2 * l * e1 * e7 - 0.72e2 * A * e2 * l * e5 * e7 - 0.72e2 * A * e6 * l * e1 * e7 + 0.72e2 * A * e6 * l * e5 * e7 + 0.18e2 * A * l ^ 2 * e8 * e2 * e4 + 0.6e1 * A * l ^ 2 * e8 * e1 * e3 - 0.6e1 * A * l ^ 2 * e8 * e3 * e5 - 0.18e2 * A * l ^ 2 * e8 * e4 * e6 + 0.6e1 * A * l ^ 2 * e4 * e1 * e7 - 0.6e1 * A * l ^ 2 * e4 * e5 * e7 + 0.6e1 * A * e2 * l ^ 2 * e3 * e7 - 0.6e1 * A * e6 * l ^ 2 * e3 * e7 + 0.24e2 * A * l ^ 3 * e8 ^ 3 - 0.3e1 * A * l ^ 3 * e4 ^ 3 - 0.42e2 * A * e2 * l ^ 2 + 0.42e2 * A * e6 * l ^ 2 - 0.324e3 * A * e2 ^ 2 * e6 + 0.108e3 * A * e2 * e1 ^ 2 + 0.108e3 * A * e2 * e5 ^ 2 + 0.324e3 * A * e2 * e6 ^ 2 - 0.108e3 * A * e6 * e1 ^ 2 - 0.108e3 * A * e6 * e5 ^ 2 + 0.14e2 * A * l ^ 3 * e4 + 0.1680e4 * II * l * e4 + 0.3360e4 * II * l * e8 - 0.216e3 * A * e2 * e1 * e5 + 0.216e3 * A * e6 * e1 * e5 + 0.6e1 * A * l ^ 3 * e8 * e4 ^ 2 + 0.2e1 * A * l ^ 3 * e8 * e3 ^ 2 - 0.9e1 * A * l ^ 2 * e8 ^ 2 * e2 - 0.9e1 * A * l ^ 3 * e8 ^ 2 * e4 + 0.9e1 * A * l ^ 2 * e8 ^ 2 * e6 + 0.36e2 * A * l * e8 * e1 ^ 2 + 0.36e2 * A * l * e8 * e5 ^ 2 + 0.108e3 * A * l * e8 * e2 ^ 2 + 0.108e3 * A * l * e8 * e6 ^ 2 - 0.3e1 * A * l ^ 3 * e4 * e7 ^ 2 + 0.9e1 * A * l ^ 2 * e4 ^ 2 * e2 - 0.3e1 * A * l ^ 3 * e4 * e3 ^ 2 - 0.9e1 * A * l ^ 2 * e4 ^ 2 * e6 - 0.3e1 * A * e2 * l ^ 2 * e7 ^ 2 + 0.3e1 * A * e2 * l ^ 2 * e3 ^ 2 + 0.3e1 * A * e6 * l ^ 2 * e7 ^ 2 - 0.3e1 * A * e6 * l ^ 2 * e3 ^ 2 + 0.24e2 * A * l ^ 3 * e8 * e7 ^ 2 + 0.108e3 * A * e2 ^ 3 - 0.108e3 * A * e6 ^ 3 - 0.5040e4 * II * e6 + 0.5040e4 * II * e2) / 0.840e3;


beta.L13   = [beta1;beta2;beta3;beta4];
beta.L23   = [beta5;beta6;beta7;beta8];