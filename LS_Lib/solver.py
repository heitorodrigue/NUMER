import math

class LSSolver:
    def __init__(self,augMat):
        self.augMat = augMat

    def print_matrix(self, matrix=None, format_spec="7.2"):
        if matrix is None:
            matrix = self.augMat

        for row in matrix:
            for value in row:
                print(f"{value:{format_spec}f}", end= " ")
            print()

    def fronebius_norm(self, matrix=None):
        """
        Calcula a norma de Frobenius da matriz fornecida.
        Se nenhuma matriz for fornecida, calcula para self.augMat.  
        
        :param matrix: matriz (lista de listas) opcional
        :return: norma de Frobenius (float)
        """
        if matrix is None:
            matrix = self.coef_matrix()

        total = 0.0
        for row in matrix:
            for value in row:
                total += value ** 2

        return math.sqrt(total)

    @staticmethod
    def zeros_matrix(rows, cols):
        """
        Cria uma matriz preenchida com zeros
        :param rows: numero de linhas
        :param cols: numero de colunas
        :return: matriz de zeros (list of lists)
        """
        return [[0.0 for _ in range(cols)] for _ in range(rows)]

    def coef_matrix(self):
        """
        Recupera a matriz de coeficientes (sem a coluna de termos independentes)
        :return: Matriz de coeficientes (list of lists)
        """
        # Recuperando as dimensoes da matriz aumentada
        rows = len(self.augMat)
        cols = len(self.augMat)

        # Criando nova matriz de zeros para os coeficientes
        MC = self.zeros_matrix(rows, cols - 1)

        # Copiando os coeficientes (tudo exceto a ultima coluna)
        for i in range(rows):
            for j in range(cols - 1):
                MC[i][j] = self.augMat[i][j]

        return MC

    def _pivot_check_and_swap(self,row_index):
        """
        Verifica se o pivo na linha 'row_index' eh zero.
        Se for, tenta trocar com uma linha abaixo que tenha pivo nao-nulo.
        Lanca um erro se a matriz for singular.
        """        
        n = len(self.augMat)

        if self.augMat[row_index][row_index] == 0:
            for j in range(row_index + 1, n):
                if self.augMat[j][row_index] != 0:
                    self.augMat[row_index], self.augMat[j] = self.augMat[j], self.augMat[row_index]
                    return # troca feita com sucesso
                raise ValueError("Matriz singular! Nao eh possivel aplicar o metodo(pivo zero).")
            
    def gaussian_elimination(self):
        """
        Executa o metodo de Eliminacao de Gauss para triangular a matriz aumentada.
        Modifica self.augMath on-place.
        """            
        n = len(self.augMat)
        m = len(self.augMat[0])
        num_swaps = 0

        for i in range(n):
            # Atencao: divide pelo pivo apenas para calcular coeficiente correto
            coef = self.augMat[j][i] / pivot
            # O laco a seguir pode otimizar comecando de i, ja que antes eh zero
            for k in range(i, m):
                self.augMat[j][k] -= coef * self.augMat[i][k]

        return num_swaps

    def determinantt(self):
        # Extrair matriz coeficiente
        coef = self.coef_matrix()

        # Salvar a matriz original
        original_augMat = self.augMat

        # Substituir por copia para nao modificar o original
        self.augMat = [row[:] for row in coef]

        # Garantir matriz quadrada
        n = len(self.augMat)
        if any(len(row)!= n for row in self.augMat):
            self.augMat = original_augMat
            raise ValueError("Matriz de coeficientes nao eh quadrada. Determinante indefinido.")

        # Executar eliminacao para triangularizacao e contar swaps
        num_swaps = self.gaussian_elimination()

        # Calcular produto da diagonal
        det = 1
        for i in range(n):
            det *= self.augMat[i][i]

        # Ajustar sinal pelo numero de swaps
        if num_swaps % 2 == 1:
            det = -det

        # Copiar a matriz triangular para retornar
        triangular = [row[:] for row in self.augMat]

        # Restaurar matriz original
        self.augMat = original_augMat

        return det, triangular

    def gauss_jordan_elimimation(self):
        """
        Executa o metodo de Gauss-Jordan para transformar a matriz aumentada
        em forma reduzida por linhas (Row Reduced Echelon Form - RREF).
        Modifica self.augMat in-place.
        """

        n = len(self.augMat)
        m = len(self.augMat[0])

        for i in range(n):
            self._pivot_check_and_swap(i)

            pivot = self.augMat[i][i]
            for k in range(m):
                self.augMat[i][k] /= pivot

            for j in range(n):
                if i != j:
                    coef = self.augMat[j][i]
                    for k in range(m):
                        self.augMat[j][k] -= coef * self.augMat[i][k]

    def inverse(self):
        """
        Calcula a inversa da matriz de coeficientes usando Gauss-Jordan.
        :return: Matriz inversa (list of lists)
        :raises ValueError: se a matriz nao for inversivel
        """
        A  = self.coef_matrix()
        n = len(A)

        #Verificar se A eh quadrada
        if any(len(row) != n for row in A):
            raise ValueError("A matriz de coeficientes nao eh quadrada. Inversa indefinida.")

        #Verificar se o determinante eh zero 
        det, _ = self.determinantt()
        if det == 0:
            raise ValueError("Matriz nao eh inversivel(det zero).") 

        # Criar identidade
        I = [[float(i == j) for j in range(n)] for i in range(n)]

        #Construir matriz aumentada [A | I]
        augmented = [A[i] + I[i] for i in range(n)]

        #Salva matriz original
        original_augMat = self.augMat
        self.augMat = [row[:] for row in augmented]

        try:
            self.gauss_jordan_elimimation()
        except ValueError:
            self.augMat = original_augMat
            raise

        #Extrair a parte direita(matriz inversa)
        inverse_mat = [row[:] for row in self.augMat]

        #Restaurar matriz original
        self.augMat = original_augMat

        return inverse_mat

    def condition_number(self):
            """
            Calcula o numero de condicao da matriz de coeficientes com base na norma de
            Frobenius e do valor da inversa da matriz.
        
            :return: numero de condicao (float)
            :raises ValueError: se a matriz nao for inversivel
            
            """
            norm_A = self.fronebius_norm() #||A||_F

            try:
                inverse_A = self.inverse() #A-1
            except ValueError as e:
                raise ValueError(f"Nao eh possivel calcular o numero de condicao: {e}")

            norm_inv_A = self.fronebius_norm(inverse_A) # ||A-1||_F

            return norm_A * norm_inv_A

    def validate_system(self, condidition_threshold=1e2):
        """
        Verifica se o sistema linear eh bem posto e estavel numericamente.
        Criterios:
        - Matriz quadrada e determinante diferente de zero (sistema bem posto)
        - Numero de condicao nao excessivamente alto (evita instabilidade)
        :param condition_threshold: limite aceitavel para o numero de condicao (
        default: 100)
        :return: dicionario com status e mensagens
        :raises ValueError: se o sistema nao for bem posto
        """
        report = {}

        # Verificar se a matriz de coeficientes eh quadrada
        A = self.coef_matrix()
        n = len(A)
        if any(len(row) != n for row in A):
            raise ValueError("Matriz de coeficientes nao eh quadrada. Sistema mal posto.")

        #Verificar determinante
        det, _ = self.determinantt()
        report["determinante"] = det

        if det == 0:
            raise ValueError("Determinante nulo. Sistema nao possui solucao unica (mal posto).")
        else:
            report["well_posed"] = True

        #Calcular numero de condicao
        kappa = self.condition_number()
        report["condition number"] = kappa

        if kappa > condidition_threshold:
            report["sensitivity_warning"] = True
            report["message"] = (
                f"Sistema eh sensivel a perturbacoes. "
                f"\nNumero de condicao alto (k = {kappa:.2f} > {condition_threshold :.0f})."
            )
        else:
            report["sensitivity_warning"] = False
            report["message"] = (
                f"Sistema considerado numericamente estavel. "
                f"k = {kappa:.2f}"
            )

        return report

            
