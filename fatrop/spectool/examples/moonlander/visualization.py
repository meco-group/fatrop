# visualize the Ftot_sol vectors using matplotlib quiver origin given by tranformation matrices F_r_sol, matrices along first axis
import matplotlib.pyplot as plt
def visualize(F_r_sol, Ftot_sol, F1_sol, F2_sol):
    plt.quiver(F_r_sol[:,0,2], F_r_sol[:,1,2], Ftot_sol[:,0], Ftot_sol[:,1])
    plt.show()

    # visualize F1_sol and F2_sol using matplotlib
    plt.plot(F1_sol)
    plt.plot(F2_sol)
    plt.show()
