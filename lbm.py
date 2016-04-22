import numpy as np


def init():
    pass


def streaming():
    pass


def collision():
    pass


def equilibrium(p, e, u, w, c):
    eu = np.dot(e, u.transpose(1, 0, 2)) / c
    uu = 3. / 2. * (u[0] ** 2 + u[1] ** 2) / c
    s = np.empty((len(e), u.shape[1], u.shape[2]))
    for i in range(len(e)):
        s[i, :, :] = w[i] * (3. * eu[i] + 9. / 2. * eu[i] ** 2 - uu)
    feq = np.empty_like(s)
    for i in range(len(e)):
        feq[i] = w[i] * p + s[i] * p
    return feq


def density(f):
    return np.sum(f, axis=0)


# TODO should be multiplied with c???
def velocity(p, f, e, c):
    return np.dot(e.transpose(), f.transpose((1, 0, 2))) / p


def lbm():
    c = 1.0
    e = np.array([(x, y) for x in [-1, 0, 1] for y in [-1, 0, 1]])
    w = np.full(9, 1. / 36.)
    for idx, ei in enumerate(e):
        norm = np.linalg.norm(ei)
        if norm == 0.:
            w[idx] = 4. / 9.
        elif norm == 1.:
            w[idx] = 1. / 9.

    f = np.random.rand(9, 5, 5)
    f[:, :, 0] = 2.
    p = density(f)
    u = velocity(p, f, e, c)
    # print("f=\n" + str(f))
    # print(e)
    # print("p =\n" + str(p))
    feq = equilibrium(p, e, u, w, c)
    # print("feq =\n" + str(feq))


def main():
    print("hi there!")
    lbm()


if __name__ == "__main__":
    main()
