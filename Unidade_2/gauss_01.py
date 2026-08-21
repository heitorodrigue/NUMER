from LS_Lib.solver import LSSolver

if __name__ == "__main__":
    matrix = [
        [1, 2, 3],
        [4, 5, 6],
        [7, 8,]
    ]

    print("Analisando o sistema linear...\n")

    lss = LSSolver(matrix)

    try:
        report = lss.validate_system()
        print("Avaliacao do sistema linear proposto.")
        print(f"Determinante: {report['determinant']:.4f}")
        print(f"Numero de condicao (k): {report['condition_number']:.2f}")
        print(f"Aviso de sensibilidade numerica: {report['sensitivity_warning']}")
        print(f"Mensagem: {report['message']}\n")

        print("Matriz aumentada origianl:")
        lss,print_matrix()

        lss.gauss_jordan_elimination

        print("\n Matriz na forma reduzida (RREF):")
        lss.print_matriz()

        solution = [row[-1] for row in lss.augMat]
        print("\n Solucao do sistema:")
        for i, val in enumerate(solution):
            print(f"x{i + 1} = {val:.4f}")

    except ValueError as e: 
        print(" Erro durante a resolucao do sistema:")
        print(str(e))
