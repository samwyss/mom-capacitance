from Solver import Solver


def main():
    solver = Solver(10, 1e-2, 1, 0.1)
    solver.solve()


# call main function after declaration
if __name__ == "__main__":
    main()
