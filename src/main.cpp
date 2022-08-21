

void Demo_Vector();
void Demo_VectorN();
void Demo_Matrix();
void Demo_MatrixNM();
void Demo_Tensors();
void Demo_Geometry();
void Demo_CoordTransf();
void Demo_CoordSystem();

void Demo_Polynom();
void Demo_Function();

void Demo_Derivation();
void Demo_Field_operations();

void Demo_Integration();
void Demo_Interpolators();
void Demo_Interpolated_Function();

void Demo_LinearAlgEqSolvers();
void Demo_DiffEqSolvers();

int main(int, char**) 
{
    // Demo_CoordSystem();
    // Demo_CoordTransf();
    // Demo_Derivation();
    // Demo_DiffEqSolvers();
    Demo_Field_operations();
    // Demo_Function();
    // Demo_Geometry();
    // Demo_Integration();
    // Demo_Interpolated_Function();
    // Demo_Interpolators();
    // Demo_LinearAlgEqSolvers();
    // Demo_Matrix();
    // Demo_MatrixNM();
    // Demo_Polynom();
    // Demo_Tensors();
    // Demo_Vector();
    // Demo_VectorN();


/* TODO
    1. VectorN osnovne operacije
    2. Matrix formatirani ispis
        - da se može inicijalizirati s vector<vector<double>>
        - inicijalizacija s initializer list?
    3. Matrix - LU rješavanje
    4. Matrix - eigenvalues & vectors
    5. Interpolacija - 
        - model TabulatedFunction
            - TabulatedFunction 2D, 3D, 4D
        - polint
        - rational polint
        - spline
    6. CoordinateSystems
        - Point3D, Point3DSpherical, Cylindrical, 
        - constexpr za izračune izraza transformacije koordinata?
        - VectorN ... sve varijante
        - METRIKA!!!

- integracija
    - multidimensional integration over regular domains
    - usporediti integraciju u Kartezijevim i sfernim koord nekog volumnog integrala preko kugle!
- root finding
- diff eq 

- i onda, rješenje diff Maxwella npr je neka vektorska funkcija
    - i onda idem izračunati Gaussa i Stokesa numerički

    - usporedba integracije diff jedn s poznatim rezultatom, 
*/

}