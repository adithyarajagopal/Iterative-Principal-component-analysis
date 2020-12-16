import untitled17
import numpy as np
from numpy.linalg import cholesky
from numpy.linalg import inv
from numpy.linalg import svd
from numpy import transpose
import math
from scipy.stats import chi2
#FUNCTION MINIMIZE_KRONECKER INPUTS:Y1 : DATA MATRIX, A1: ESTIMATED CONSTRAINT MATRIX, ROWIND: ROWINDEX(ROW INDEX OF NONZERO ERROR COVARIANCE MATRIX ELEMENTS),COLIND: COLUMN INDEX(COLUMN INDEX OF NON ZERO ERROR COVARIANCE MATRIX ELEMENTS)
#                            OUTPUTS: SIGMA(ERROR COVARIANCE MATRIX)
#THE SIZE OF EITHER ROW INDEX OR COLUMN INDEX CANNOT EXCEED m(m+1)/2 WHERE m IS THE NUMBER OF CONSTRAINTTS IN THE MODEL
#FUNCTION CALLS(DEPENDENCIES):VARIANCE_MAPPER_ARRAY_TO_SIGMA:INPUTS:ARRAY(TRANPOSE OF VEC(ERROR COVARIANCE MATRIX), ROWINDEX, COLUMN INDEX)
# 
                                        #   OUTPUTS:SIGMA(ERROR COVARIANCE MATRIX)
def minimize_kronecker(Y1, A1, rowind1, colind1):
    import numpy as np
    from numpy.linalg import inv
    
    COUNT_SUM = 0
    (n,N) = Y1.shape
    COUNT_MATRIX = []
    COUNT_MATRIX = np.array(COUNT_MATRIX, dtype = float)
    
    #COMPUTING ESTIMATE OF COVARIANCE OF RESIDUALS USING ESTIMATED CONSTRAINT MATRIX
    #COUNT_SUM IS THE MATRIX OF COVARIANCE OF RESIDUALS 
    for i in range(N):
        COUNT_MATRIX = np.dot(A1,Y1[ : ,i])
        
        COUNT_SUM = COUNT_SUM +( np.dot(COUNT_MATRIX,transpose(COUNT_MATRIX))/N)
   
    LIST = []
    r,c = A1.shape
    k = n
    k1 = len(rowind1)
    for i in range(k1):
        for j in range(k1):
            if(i == j):
               
                K = k*(colind1[i]) + (rowind1[i])
                LIST.append(K)
    #THE FOLLOWING INBUILT NUMPY FUNCTION COMPUTES THE KRONECKER PRODUCT BETWEEN TWO MATRICES, 
    #HERE G IS AN (M**2,N**2) MATRIX WHERE M IS THE NUMBER OF CONSTRAINTS AND N IS THE NUMBER OF PROCESS VARIABLES
    G = np.kron(A1,A1)
   
    #LIST IS THE COLUMNS CORRESPONDING TO THE NON ZERO ELEMENTS OF THE ERROR COVARIANCE MATRIX
    G = G[ : ,LIST]
   
    
    resid = COUNT_SUM
   
    vec_resid = transpose(transpose(resid).flatten())
    
    
    
    
    vec_sigma = np.dot(inv(np.dot(transpose(G),G)),np.dot(transpose(G),vec_resid))
    
    vec_sigma_t = transpose(vec_sigma)
    vec_sigma_t = np.array(vec_sigma_t, dtype = float)
    vec_sigma_t = vec_sigma_t.reshape(1,len(rowind1))
   
    vec_sigma_t = np.absolute(vec_sigma_t)
    print(vec_sigma_t)
    
    sigma1 = Variance_mapper_array_to_sigma(vec_sigma_t, rowind1, colind1)
    sigma1 = sigma1.reshape(n,n)
    print(sigma1.shape)
   
    
   
 
    return sigma1
#THIS FUNCTION IS NOT REQUIRED            
def Variance_Mapper_sigma_to_array(sigma,rowind,colind):
    rowind = np.array(rowind, dtype = float)
    colind = np.array(colind, dtype = float)
    rowind = rowind.sorted()
    colind = colind.sorted()
    list_ = [] 
    k = len(rowind)
    for i in range(k):
        for j in range(k):
            if(i == j):
                 list_.append(sigma[rowind[i], colind[j]])
            
    print(list_)       
    return list_
#THIS FUNCTION CONVERTS AN ARRAY TO A MATRIX BY MAPPING THE ELEMENTS CORRESPONDING TO THE ROW AND COLUMN INDEX
#FUNCTION INPUTS AND OUTPUTS EXPLAINED ABOVE MINIMIZE_KRONECKER FUNCTION   
def Variance_mapper_array_to_sigma(array1, rowind, colind):
    rowind = np.array(rowind, dtype = int)
    colind = np.array(colind, dtype = int)
    n = 5
    
    #CREATING AN N*N MATRIX OF ZEROES
    sigma = [0 for i in range(n) for j in range(n)]
    sigma = np.asmatrix(sigma, dtype = float)
    sigma = sigma.reshape(n,n)
    array1 = np.array(array1, dtype = float)
    k = len(rowind)
    array1 = array1.reshape(1,k)
    #INSERTING THE ELEMENTS OF THE ARRAY INTO THE MATRIX AND RETURNS THE RESULTING MATRIX
    for i in range(k):
        for j in range(k):
            if(i == j):
               
                rowind[i] = int(rowind[i])
                colind[j] = int(colind[j])
               
               
                sigma[rowind[i], colind[j]] = array1[0,i]
                
    
    print(sigma)
    return sigma
# SVD_MATRIX FUNCTION INPUTS: DATAMATRIX, NUMBER OF CONSTRAINTS
#                     OUTPUTS: SUM(SUM OF THE LAST THREE SINGULAR VALUES), ESTIMATED CONSTRAINT MATRIX
    
def svd_matrix(y_matrix, m):
   
    U,S,VT = svd(y_matrix)
    (n,N) = y_matrix.shape
    SUM = S[n - m:]
    print("ARRAY OF LAST THREE SINGULAR VALUES 12345")
    print(SUM)
    #THIS IS BECAUSE Y = MATH.SQRT(N)(U*S*TRANSPOSE(V))
    SUM = sum(SUM)/float(math.sqrt(N))
    
    
    
    (n, N) = y_matrix.shape
    if(N>n):
        
        
        
        U1 = np.asmatrix(U)
   #     S1 = np.asmatrix(S)
   #     VT1 = np.asmatrix(VT)
        
        U1 = U1[ : , n - m: ]  
        A = transpose(U1)
        
        return(A,SUM)
#CHOLESKY-DECOMPOSITION INPUTS: ERROR COVARIANCE MATRIX
#                       OUTPUTS:INVERSE OF THE CHOLESKY MATRIX        
def Cholesky_Decomposition(sigma1):
    from numpy.linalg import cholesky
    from numpy.linalg import inv
    
    
    L = cholesky(sigma1)
  
    L1 =   inv(L)
    return(L1)
#THE FOLLOWING FUNCTION COMPUTES THE NOISY DATA MATRIX OBTAINED FROM A CSV FILE 
    #THE DATAMATRIX USED CAN BE EITHER NOISED OR DENOISED 


#XPCA FUNCTION INPUTS:SIGMA:ERROR COVARIANCE MATRIX,ROWIND: ROWINDEX, COLIND: COLUMN INDEX, M : NUMBER OF CONSTRAINTS, TYPE OF MODEL, TYPE OF OPTIMISER
     #         OUTPUTS:A_NEXT : ESTIMATED CONSTRAINT MATRIX
     #                 SUM_SINGULAR_VALUES: SUM OF LAST THREE SINGULAR VALUES
     #                 U_ : LEFT SINGULAR VECTOR MATRIX
     #FUNCTION DEPENDENCIES(FUNCTION CALLS):CHOLESKY_DECOMPOSITION
                                            #MINIMIZE_KRONECKER
                                            #SVD_MATRIX
                                            #THE DESCRIPTION OF THE ABOVE FUNCTIONS HAVE BEEN DONE ABOVE
                                            #THE FOLLOWING FUNCTION COMPUTES PCA, MLPCA, IPCA DEPENDING ON THE ARGUMMENT PROVIDED BY THE USER
                                            
                                            
def XPCA(DATAMATRIX,SIGMA,m_cons,row_index,col_index, type_of_model, type_of_optimiser):
    #STORES A COPY OF THE ORIGINAL MATRIX IN Y
    Y = DATAMATRIX.copy()

    if(type_of_model != "IPCA"):
        L = cholesky(SIGMA)
        (n,N) = DATAMATRIX.shape
        DATAMATRIX = np.dot(inv(L), DATAMATRIX)
        U,S,VT = svd(DATAMATRIX)
        print(S)
        print("ARRAY OF SINGULAR VALUES")
        U = U[ : ,n - m_cons: ]
        A_next = transpose(U)
        
    else:
        print("TYPE OF MODEL IS IPCA")
        (n,N) = DATAMATRIX.shape
        sigma = np.eye(n)
        LIST1_ = []
        #THE CONTROL COMES OUT OF THE LOOP IF EITHER THE BREAKING POINT IS MET OF THE NUMBER OF ITERATIONS EXCEEDS A CERTAIN VALUE
        for i in range(100):
            (n,N) = DATAMATRIX.shape
            
         #IF THE NUMBER OF CONSTRAINTS IS GIVEN THE THE FOLLOWING CODE WILL BE EXECUTED, ELSE IT COMPUTES THE NUMBER OF CONSTRAINTS USING A HEURISTIC  
            if(m_cons!= None):
           
                Linv = Cholesky_Decomposition(sigma)
                DATAMATRIX = np.dot(Linv,Y)
                A_next,SUM_SINGULAR_VALUES = svd_matrix(DATAMATRIX,m_cons)
                #APPENDS THE SUM OF LAST THREE SINGULAR VALUES IN A LIST 
                LIST1_.append(SUM_SINGULAR_VALUES)
                
                A_next = np.dot(A_next,Linv)
              
                print(LIST1_)
                
               
                if(type_of_optimiser =="LEAST SQUARES"):
            
                    sigma = minimize_kronecker(Y,A_next,row_index,col_index)
               
                else:
                    print("ERROR: ENTER A VALID OPTIMISER")
                    break
               #BREAKING POINT OF THE LOOP 
                if((i > 4)):
                    print("SUM OF LAST THREE SINGULAR VALUES")
                    U_,S_,VT_ = svd(DATAMATRIX)
                    print(U_)
                    
                    print(SUM_SINGULAR_VALUES)
                    print(A_next)
                    break
            else:
                #IF THE NUMBER OF CONSTRAINTS ARE NOT GIVEN IT IS ESTIMATED USING A HEURISTIC AND THE FUNCTION IS CALLED AGAIN
                #HERE P IS THE TOTAL NUMBER OF EIGEN VALUES
                #Q IS THE NUMBER OF PRINCIPAL COMPONENTS
                #P - Q IS THE NUMBER OF CONSTRAINTS
                 U,S,VT = svd(Y)
                 print(S)
                 l = S**2
                 n1 = N#THERE MIGHT BE A PROBLEM HERE
                 p = n
                 n2 = n1 - (2*p+11)/6#IF NOT CHECK WITH MATLAB CODE
                 for q in range(p - 2, 0, -1):
                     d = p - q
                     d1 = (d*(d+1))/2
                     
                     if(d1 > (n - 1)):
                         l_avg = sum(l[q + 1 : p ])/(p - q)
                         print(l_avg)
                         st = n2*((p - q)*np.log(l_avg) - sum(np.log(l[q + 1 : p ])))#THERE MIGHT BE A PROBLEM IN LOOOPING BACKWARDS
                         #df is the degrees of freedom
                         df = ((p - q + 2)*(p - q - 1))/2
                         #5 percent level of significance assumed
                         chi_ = chi2.ppf( 0.95, df)
                        
                         print(st)
                         print(chi_)
                         #IF THE TEST STATISTIC IS LESSER THAN THE CRITICAL VALUE THE NULL HYPOTHESIS CANNOT BE REJECTED
                         if(st < chi_):
                             d = p - q
                             break
                 d = p - q    
                 m_cons = d
                 print("NUMBER OF CONSTRAINTS")
                 print(d)
                
                     
          #       eig = S**2/float(100)
          #       print("ARRAY OF EIGEN VALUES")
           #      print(eig)
            #     eig = np.array(eig, dtype = float)
            #     S = sum(eig)
            #     print("SUM OF EIGEN VALUES")
            #     print(S)
                 #A.CUMSUM() CALCULATES THE CUMULATIVE SUM OF A
                 #FOR EG:
                 #      A = [1,2,3,4,5]
                 #      A.CUMSUM() = [1,3,6,10,15]
             #    eig_cumsum = eig.cumsum()
             #    print(eig_cumsum)
              #   print("CUMULATIVE SUM")
               #  for i in range(len(eig_cumsum)):
                #     percent_eig = (eig_cumsum/S)*100
                 #    print("PERCENT EIGEN VALUES")
                     
                  #   print(percent_eig)
       #              if(percent_eig[i]>float(99.99)):
        #                 print(i)
         #                break
          #       eig_new = eig[ : i+1]
           #      m = (len(eig) - i)
           #      print("M")
            #     print(m)
             #    m_cons = m
             
                 XPCA(DATAMATRIX,SIGMA,m_cons,row_index,col_index, "IPCA","LEAST SQUARES")
    return A_next,SUM_SINGULAR_VALUES,U_  
#TESTING_FUNCTION INPUTS:A_NEXT :ESTIMATED CONSTRAINT MATRIX,A_TRRUE: TRUE CONSTRAINT MATRIX OBTAINED FROM FIRST PRINCIPLE EQUATIONS
#                 OUTPUTS: ALPHA      
def testing_function(A_true,A_next):
    from numpy.linalg import inv
    from numpy import linalg as LA
    import numpy as np
   
    alpha = []
   
    print(A_next.shape)
    
    print(type(alpha))
    m_cons = A_next.shape[0]
    for i in range(m_cons):
        
        term1 = A_true[i, : ]
        print(term1.shape)
        term1 = np.array(term1, dtype = float)
        
        term2 = np.dot(A_true[i,:],np.dot(A_next.T,np.dot(inv(np.dot(A_next,A_next.T)),A_next)))
        term2 = np.array(term2, dtype = float)
       
        TERM = term1 - term2
        print(TERM)
     
        #IF ANY VECTOR IN A_NEXT LIES IN THE HYPERPLANE FORMED BY A_TRUE, THEN IT WILL NOT BE APPENDED TO THE LIST
        if((term1 - term2).all() != None):
         
        
            alpha.append(TERM)
    
    print(alpha)
  
    for j in range(len(alpha)):
        #COMPUTES THE MAGNITUDE OF THE LIST
        
        alpha[j] = LA.norm(alpha[j])
    print("ARRAY OF ALPHA I'S")    
    print(alpha)
    print('ALPHA')
    SUM = sum(alpha)
    
    return SUM
def main():
    
    m_cons1 = None
    
    print(m_cons1)
    
    row_index1 = [0,1,2,3,4]
    col_index1 = [0,1,2,3,4]
    SIGMA1 = np.eye(5)#FOR PCA AND IPCA,
    
    A_true1 = np.matrix([[1,1,-1,0,0],[0,0,1,-1,0],[0,-1,0,1,-1]])
    DATAMATRIX1 = untitled17.DATAMATRIX
    A_next1,SUM_,U1_ = XPCA(DATAMATRIX1,SIGMA1,m_cons1,row_index1,col_index1,"IPCA", "LEAST SQUARES")
    #ONLY IF THE ESTIMATED CONSTRAINT MATRIX AND THE TRUE CONSTRAINT MATRIX HAVE THE SAME SHAPE THE SHORTEST DISTANCE BETWEEN EACH VECTOR AND THE SUBSPACE FORMED BY A_TRUE CAN BE FOUND
    if(A_true1.shape == A_next1.shape):
        alpha = testing_function(A_true1,A_next1)
        print(alpha)
    print("CONVERGED SINGULAR VALUE SUM")    
    print(SUM_) 
    
    
if __name__ == "__main__":
    main()
         
    
    

    
                     
                 
