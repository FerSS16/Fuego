/* ------------------------------ INCLUYE ----------------------------------- */
 //#include<stdio.h>
 #include <cstdlib>
 #include <math.h>
 #include <string.h>
 //#include<conio.h>
 //#include<dos.h>
 #include <time.h>

 #include <iostream>
 #include <ctime>
 #include <cstdlib>
 #include <iomanip>
 #include <fstream>
 #include <string>
 #include <vector>
 #include <map>

 using namespace std;

/*----------------------------------------------------------------------------*/
///cantidad de filas y columnas de la matriz
 const int fil = 100;
 const int col = fil;

///cantidad de veces que se recorre la matriz de estados
 const int pasos = 100;
 const int salteo = 0;
 const int repeticiones = 1;

///variables
 //std::vector<float> valores_de_P;// PROBABILIDAD DE CRECIMIENTO
 //double valores_de_P;
  const float P=0.1;
 
 const float g = 0;// INMUNIDAD
// const float h = P/10;// FRECUENCIA DE GENERACION ESPONTANEA DE FUEGOS (F)

///raster scan (H-K)
#define min(a,b) (a>b?b:a)
#define max(a,b) (a>b?a:b)

///datos del conteo de cajas (dimension fractal)
#define NCAJAS 32 /// n�mero de cajas para el conteo
#define LMIN 1 /// tama�o m�nimo de la caja
#define LMAX 64 /// tama�o m�ximo de la caja

/// n�mero de cajas que contienen sitios con fuego
int ncajas[NCAJAS];
//int i, j, k;
double L, N, logL, logN, D;
//FILE *fp;

///objetos que se van a usar
 //matrices de probabilidades
 float Matriz_probabilidades[fil][col];
 float Matriz_inmunity[fil][col];
 float Matriz_f[fil][col];

 //Evolucion
 float Matriz_estados[fil][col];
 float Matriz_pasos[fil][col];
 float Matriz_ayuda[fil][col];

 char Matriz_moleculas[fil][col];
 char H = H;
 char C = C;
 char O = O;

 //Hoshen-Kopelman
 bool Matriz_occupied[fil][col];
 int label[fil][col];
 const int MAX_SIZE = fil;
 int largest_label;
 vector<int> labels(MAX_SIZE * MAX_SIZE);

 //condiciones periodicas
 int sum[fil];
 int res[fil];

 //cantidades
 double fuegos;
 double vacios;
 double arboles;
 double fvt;
 double vvt;
 double avt;

 //dimension fractal
 double radio;

 ///variable entera que indica un radio entre un valor i+1 e i
 int var;
 const int NHS=fil;    // N�meros de puntos en el Histograma
 int Histograma[NHS+1];
 int mediciones=0;
 int a;
 int b;

///*----------------------------------------------------------------------------*/
///*----------------------GENERADOR RANDOM--------------------------------------*/
///*----------------------------------------------------------------------------*/

 const double long_max=2147483648.0;
 const double dos_long_max=4294967296.0;
 unsigned int condi_55(int);
 int mas1_55[55],mas31_55[55];
 long int randi[55];
 int j_ran,j_ran2; ///las usa el generador de random

/*----------------------------------------------------------------------------*/
void randomizar(void)
{   int i_ran,j_ran2;
    j_ran=0;
    for(i_ran=0;i_ran<55;i_ran++) randi[i_ran]=rand();
    for(i_ran=0;i_ran<1000;i_ran++)
      for(j_ran2=0;j_ran2<1000;j_ran2++)
       {j_ran=mas1_55[j_ran];
randi[j_ran]=randi[j_ran]+randi[mas31_55[j_ran]]; }

}
/*----------------------------------------------------------------------------*/

double ran01(void)
{   j_ran=mas1_55[j_ran];
    randi[j_ran]=randi[j_ran]+randi[mas31_55[j_ran]];
    return (((double)randi[j_ran]+long_max)/dos_long_max);
}

/*----------------------------------------------------------------------------*/

unsigned int condi_55(int z)
{   int w;
    if ((z>=0)*(z<55))  w=z;
    if (z<0)            w=55+z;
    if (z>=55)          w=z-55;
    return w;
}

/*----------------------------------------------------------------------------*/

void contorno(void)
{for(j_ran2=0 ; j_ran2<55 ; j_ran2++)
   { mas1_55[j_ran2]=condi_55(j_ran2+1);
    mas31_55[j_ran2]=condi_55(j_ran2+31); }
}

///*----------------------------------------------------------------------------*/
/*-----------------------------CONDICIONES DE CONTORNO------------------------*/

void Condiciones_Periodicas(void){

 int z;

       for(z=0;z<fil;z++){

        sum[z]=z+1;

            }

    for(z=0;z<fil;z++){

        res[z]=z-1;

            }

        res[0]=fil-1;
        sum[fil-1]=0;
 }

///*----------------------------------------------------------------------------*/
///*--------------------------MATRICES DE VARIABLES-----------------------------*/
///*----------------------------------------------------------------------------*/

void Matriz_Probabilidades(void){

     for(int f=0;f<fil;f++){
            for(int c=0;c<col;c++){
                    Matriz_probabilidades[f][c] = (float) ran01();
            }
      }
}

/*----------------------------------------------------------------------------*/

void Matriz_Inmunity(void){

     for(int f=0;f<fil;f++){
            for(int c=0;c<col;c++){
                    Matriz_inmunity[f][c] = (float) ran01();
            }
      }
}

/*----------------------------------------------------------------------------*/

void Matriz_F(void){

     for(int f=0;f<fil;f++){
            for(int c=0;c<col;c++){
                    Matriz_f[f][c] = (float) ran01();
            }
      }
}

/*----------------------------------------------------------------------------*/

void Matriz_Moleculas(void){

     for(int f=0;f<fil;f++){
            for(int c=0;c<col;c++){
                Matriz_moleculas[f][c] = 'x';
                }
        }
}

/*----------------------------------------------------------------------------*/
/*----------------------------MATRICES DE EVOLUCION---------------------------*/
/*----------------------------------------------------------------------------*/

void Matriz_Estados(void){

     for(int f=0;f<fil;f++){
            for(int c=0;c<col;c++){
               Matriz_ayuda[f][c]= (float) ran01();///matriz inicial
        }
}

///inicia con 15% de fuego y 25% de arboles

        for(int f=0;f<fil;f++){
            for(int c=0;c<col;c++){
                if(Matriz_ayuda[f][c]>0.75){Matriz_estados[f][c] = 1;}
                    else if(Matriz_estados[f][c]>0.15 && Matriz_estados[f][c]<=0.75){Matriz_estados[f][c]=0;}
                        else if(Matriz_ayuda[f][c]<=0.15){Matriz_estados[f][c]=2;}
        }
    }

///imprime matriz estados
/*
        for(int f=0;f<fil;f++){
            for(int c=0;c<col;c++){
                cout<<Matriz_estados[f][c]<<"  ";
            }
                cout<<endl;
            }
                cout<<endl;
*/
}

/*--------------------------------------------------------------*/

void Matriz_Pasos(void){

     for(int f=0;f<fil;f++){
            for(int c=0;c<col;c++){
               Matriz_pasos[f][c]=1;//matriz inicial
                }
        }
}

/*--------------------------------------------------------------*/
/*-----------------CONTADOR DE F,V Y A--------------------------*/
/*--------------------------------------------------------------*/

void Distribucion(void){

        for(int f=0;f<fil;f++){
       for(int c=0;c<col;c++){

           if(Matriz_estados[f][c]==2){
///conteo de fuegos
                fuegos+=1;
                fvt+=1;
                
///histograma de fuegos
/*
           for(int m=0;m<fil;m++){
                    for(int n=0;n<col;n++){
                        if(Matriz_estados[m][n]==2){
                            if(fuegos>0){
                                if(sqrt((f-m)*(f-m))<=(fil/2)){a=sqrt((f-m)*(f-m));}
                                if(sqrt((f-m)*(f-m))>(fil/2)){a=fil-sqrt((f-m)*(f-m));}
                                if(sqrt((c-n)*(c-n))<=(fil/2)){b=sqrt((c-n)*(c-n));}
                                if(sqrt((c-n)*(c-n))>(fil/2)){b=fil-sqrt((c-n)*(c-n));}
                            }
                                radio = sqrt((a*a)+(b*b));
                                var = ((NHS*radio)/fil);
                                Histograma[var]+=1;
                         }
                    }
               }
*/
            }
            else if(Matriz_estados[f][c]==0){
                vacios+=1;
                vvt+=1;
                
                }
            else if(Matriz_estados[f][c]==1){
                arboles+=1;
                avt+=1;
                }
         }
    }
}

/*--------------------------------------------------------------*/
/*-----------------DIMENSION FRACTAL (no me acuerdo como hace)--*/
/*--------------------------------------------------------------*/

void DF(void){

/// inicializar el arreglo de n�mero de cajas con ceros
//cout<<"log10(L)       log10(N)"<<endl;

/// calcular el n�mero de cajas que contienen sitios con fuego para diferentes tama�os de cajas
/*    for (int k = 0; k < NCAJAS; k++) {
        L = LMIN + k * (LMAX - LMIN) / (NCAJAS - 1); // tama�o de la caja
        for (int i = 0; i < L; i ++){
            for (int j = 0; j < L; j ++) {
                if (Matriz_estados[i][j] == 2) {
                    ncajas[k]++;
            }
        }
    }
        cout<<L<<"       "<<ncajas[k]<<endl;
}
*/

    for(int k = 0; k < NCAJAS; k++) {
/// tama�o de la caja
    L = LMIN + k * (LMAX - LMIN) / (NCAJAS - 1);
        for(int i = 0; i < fil; i +=L){
            for(int j = 0; j < fil; j +=L){
                if (Matriz_estados[i][j] == 2){
                    ncajas[k]++;
            }
        }
    }
//cout<<L<<"       "<<ncajas[k]<<endl;
}
/*
D = (log10(ncajas[NCAJAS-1]) - log10(ncajas[0])) / (log10(LMAX) - log10(LMIN));
cout<<ncajas[NCAJAS-1]<<endl;

/// calcular la dimensi�n fractal a partir de los datos de conteo de cajas
/// archivo para guardar los resultados
 fp = fopen("fractal.txt", "w");
 fprintf(fp, "log(L)\tlog(N)\n");
 cout<<"log(L)       log(N)"<<endl;

    for (int k = 0; k < NCAJAS; k++) {
/// tama�o de la caja
        L = LMIN + k * (LMAX - LMIN) / (NCAJAS - 1);
        logL = log10(L);
        logN = log10(ncajas[k]);
         fprintf(fp, "%.4f\t%.4f\n", logL, logN);
         cout<<logL<<"       "<<logN<<endl;
         cout<<L<<"       "<<ncajas[k]<<endl;
}
     fclose(fp);
     D = (log10(ncajas[NCAJAS-1]) - log10(ncajas[0])) / (log10(LMAX) - log10(LMIN)); // dimensi�n fractal: pendiente log log
     printf("Dimensi�n fractal: %.4f\n", D);
     cout<<endl<<"DF = "<<D<<endl;
*/
}

/*--------------------------------------------------------------*/

void Matriz_Occupied(void){

     for(int f=0;f<fil;f++){
        for(int c=0;c<col;c++){
               Matriz_occupied[f][c]=0;
          }
     }
}

/*--------------------------------------------------------------*/

void Labels(void){
     for(int x=0;x<fil*col;x++){
               labels[x]=x;
     }
}

/*--------------------------------------------------------------*/

void Label(void){

     for(int f=0;f<fil;f++){
            for(int c=0;c<col;c++){
                    label[f][c] = 0;
            }
        }
}

/*--------------------------------------------------------------*/
/*-------------------------HOSHEN-KOPELMAN----------------------*/
/*--------------------------------------------------------------*/

///FUNCIONES QUE USAR� RASTER SCAN

///encuentra marcas
int find_label(int x) {
    int y = x;

    while (labels[y] != y)
        y = labels[y];

    while (labels[x] != x)  {
        int z = labels[x];
        labels[x] = y;
        x = z;
    }

    return y;
}

///une marcas
void union_labels(int x, int y) {
    labels[find_label(x)] = find_label(y);
}

///RASTER SCAN
void raster_scan(int col, int fil, bool Matriz_occupied[MAX_SIZE][MAX_SIZE]) {

///empiezo por el primer lugar [0][0]

         if (Matriz_occupied[0][0]) {

                if (label[0][0]== 0) {
                    largest_label++;
                    label[0][0] = largest_label;
                    labels[largest_label] = largest_label;

                } else if (label[0][0] != 0) {
                    label[0][0] = 1;
                }
            }

///recorro solo una fila y una columna

    for (int x = 0; x < 1; x++) {
        for (int y = 1; y < fil; y++) {
            if (Matriz_occupied[0][y]) {

                int left = label[x-1][y];
                int above = label[x][y-1];

                if (above == 0) {
                    largest_label++;
                    label[0][y] = largest_label;
                    labels[largest_label] = largest_label;

                } else if (above != 0) {
                    label[0][y] = find_label(above);
                } else {
                    union_labels(max(left, above),min(left, above));
                    label[x][y] = find_label(min(left, above));
                }
            }
        }
    }

    for (int x = 1; x < col; x++) {
        for (int y = 0; y < 1; y++) {
            if (Matriz_occupied[x][y]) {

                int left = label[x-1][y];
                int above = label[x][y-1];

                if (left == 0) {
                    largest_label++;
                    label[x][y] = largest_label;
                    labels[largest_label] = largest_label;
                } else if (left != 0) {

                    label[x][y] = find_label(left);

                } else {
                     union_labels(max(left, above),min(left, above));
                    label[x][y] = find_label(min(left, above));
                }
            }
        }
    }

///recorro toda la matriz
    for (int x = 1; x < col; x++) {
        for (int y = 1; y < fil; y++) {
            if (Matriz_occupied[x][y]!=0) {

                int left = label[x-1][y];
                int above = label[x][y-1];

                if (left == 0 && above == 0) {
                    largest_label++;
                    label[x][y] = largest_label;
                    labels[largest_label] = largest_label;
                } else if (left != 0 && above == 0) {

                    label[x][y] = find_label(left);
                } else if (left == 0 && above != 0) {
                    label[x][y] = find_label(above);
                } else if(left != 0 && above != 0){
                    union_labels(max(left, above),min(left, above));
                    label[x][y] = find_label(min(left, above));
                }
            }
            else if(Matriz_occupied[x][y]==0){continue;}
        }
     }

  for(int l=0;l<fil*fil;l++){
    int w = find_label(l);
       for (int y = 0; y < fil; y++) {
        for (int x = 0; x < col; x++) {
            if (label[x][y]==l){

                label[x][y]=w;
                }
            }
        }
    }
}

//--------------------------------------------------------------------//
///medir el tama�o de cada cluster (no anda en el cluster, esta remediado mas adelante)

/*void Tamano_Cluster(){
    map<int, int> cluster;
    map<int, int> cantidad;

    for(int x=0; x<fil;x++){
        for(int y=0; y<col;y++){
            int labels = label[x][y];
            if(labels != 0){
                cluster[labels]++;
            }
        }
    }
///tabla que dice cuantos hay de cada tama�o -sin ordnar-

    for(auto& entry : cluster){
        cout<<entry.second<<endl;
    }cout<<endl;


    for(auto& entry : cluster){
        int tamano = entry.second;
        cantidad[tamano]++;
    }

///tabla que dice cuantos clusters de cada tama�o hay -ordenados-
//cout<<"Tama�o    Cantidad"<<endl;

ofstream clusters("tama�o", ios::app);
//ofstream fout("filename.txt", ios::app);
    for (const auto& entry : cantidad){
        int tamanio = entry.first;
        int cantidad = entry.second;

///guarda archivo tama�o
clusters<<"Tama�o    Cantidad"<<endl<<tamanio<< "         " <<cantidad<< endl;

        //cout<<tamanio <<"           "<<cantidad<<endl;
    }
    //cout<<endl;
clusters.close();

}*/
void Tamano_Cluster() {
    map<int, int> cluster;
    map<int, int> cantidad;

    // Asumimos que 'fil' y 'col' están definidas en otro lugar del código
    for(int x = 0; x < fil; x++) {
        for(int y = 0; y < col; y++) {
            int labels = label[x][y];
            if(labels != 0) {
                cluster[labels]++;
            }
        }
    }

    // Tabla que dice cuántos hay de cada tamaño - sin ordenar
    for(map<int, int>::iterator it = cluster.begin(); it != cluster.end(); ++it) {
        printf("%d\n", it->second);
    }
    printf("\n");

    // Tabla que dice cuántos clusters de cada tamaño hay - ordenados
    for(map<int, int>::iterator it = cluster.begin(); it != cluster.end(); ++it) {
        int tamano = it->second;
        cantidad[tamano]++;
    }

    // Guardar archivo tamaño
    FILE *clusters = fopen("tamaño.txt", "a");
    if (clusters == NULL) {
        perror("No se pudo abrir el archivo");
        return;
    }
    
    fprintf(clusters, "Tamaño    Cantidad\n");
    for(map<int, int>::iterator it = cantidad.begin(); it != cantidad.end(); ++it) {
        int tamanio = it->first;
        int cantidad = it->second;

        fprintf(clusters, "%d         %d\n", tamanio, cantidad);
    }
    
    fclose(clusters);
}

//--------------------------------------------------------------------//
///------------------HOSHEN-KOPELMAN---------------------------------///
//--------------------------------------------------------------------//

int Hoshen_Kopelman(){

Labels(); Label();


///pasa la matriz de 0, 1 y 2 a solo 0 y 1
   for(int f=0;f<fil;f++){
        for(int c=0;c<col;c++){
//ver si quiero arboles o fuegos:
            if(Matriz_estados[f][c] == 0 || Matriz_estados[f][c] == 2){
                Matriz_occupied[f][c] = 0;}
            else if(Matriz_estados[f][c] == 1){
                Matriz_occupied[f][c] = 1;

            }
        }
    }

    raster_scan(col, fil, Matriz_occupied);

/// Imprimir las matrices usadas
/*
    for (int y = 0; y < fil; y++) {
        for (int x = 0; x < col; x++) {
            cout << Matriz_occupied[x][y] << "        ";
        }
            cout << endl<<endl;
        }
cout<<endl;
    for (int y = 0; y < fil; y++) {
        for (int x = 0; x < col; x++) {
            cout << label[x][y] << "        ";
        }
            cout << endl<<endl;
        }
*/
    Tamano_Cluster();
    return 0;
}

/*--------------------------------------------------------------*/
/*--------------------EVOLUCION---------------------------------*/
/*--------------------------------------------------------------*/
void Pasos(float valor){

///abre el archivo video
//ofstream matriz_moleculas("Matriz_moleculas");
///abre archivo f_vs_t
ofstream f_vs_t("f_vs_t");
         f_vs_t<<"t      fuego   arboles    vacios"<<endl;


       for(int i=1;i<pasos+1;i++){

        //fuegos=0;

///calcula la dimension fractal
//    DF();

///Hoshen-Kopelman
if(i==pasos-1){
    Hoshen_Kopelman();
//    Tamano_Cluster();
//Distribucion();

  ///guarda hoshen kopelman:
/*
ofstream matriz_hk("Matriz_hk");
        matriz_hk<<"ITEM: TIMESTEP"<<endl<<i<<endl<<"ITEM: NUMBER OF ATOMS"<<endl<<fil*col<<".0"<<endl<<"ITEM: BOX BOUNDS pp pp"<<endl<<"0.0"<<" "<<fil-1<<".0"<<endl<<"0.0"<<" "<<col-1<<endl<<"0.0"<<" "<<"0.0"<<endl<<"ITEM: ATOMS id type x y z"<<endl;

        int a=0;
            for(int f=0;f<fil;f++){
                for(int c=0;c<col;c++){
                    ///usar la matriza label para clusters y la Matriz_moleculas para los fuegos/arboles.
                        matriz_hk<<a<<"    "<<label[f][c]<<"    "<<f<<"    "<<c<<"    "<<0<<endl;
                        a++;
                        }
                    }
matriz_hk.close();
*/
///archivo guardado

}

 fvt =0;
 avt=0;
 vvt=0;
 Distribucion();

///guarda distribucin temporal  f_vs_t
 f_vs_t<<i<<"   "<<fvt<<"   "<<avt<<"   "<<vvt<<endl;

   if(i>1 && fuegos==0){
       // cout<<"stop en i= "<<i<<endl;
        //Hoshen_Kopelman();
   break;
   }


///guarda el archivo de fuegos, vacios y arboles
///Primero le pongo nombre de atomos a los lugares

          for(int f=0;f<fil;f++){
            for(int c=0;c<col;c++){

            int matriz_e=Matriz_estados[f][c];

            switch(matriz_e){

        case 0: {Matriz_moleculas[f][c]= 'O';}break;
        case 1: {Matriz_moleculas[f][c]= 'H';}break;
        case 2: {Matriz_moleculas[f][c]= 'C';}break;

                            }
                        }
                    }
/*
        matriz_moleculas<<"ITEM: TIMESTEP"<<endl<<i<<endl<<"ITEM: NUMBER OF ATOMS"<<endl<<fil*col<<".0"<<endl<<"ITEM: BOX BOUNDS pp pp"<<endl<<"0.0"<<" "<<fil-1<<".0"<<endl<<"0.0"<<" "<<col-1<<endl<<"0.0"<<" "<<"0.0"<<endl<<"ITEM: ATOMS id type x y z"<<endl;

        int k=0;
            for(int f=0;f<fil;f++){
                for(int c=0;c<col;c++){
                    ///usar la matriza label para clusters y la Matriz_moleculas para los fuegos.
                        matriz_moleculas<<k<<"    "<<Matriz_moleculas[f][c]<<"    "<<f<<"    "<<c<<"    "<<0<<endl;
                        k++;
                        }
                    }
*/

///archivo guardado

///sortea las variables

    Matriz_Probabilidades(); Matriz_Inmunity();//Matriz_F();

///Evolucion

        for(int f=0;f<fil;f++){
                for(int c=0;c<col;c++){

        int matriz_e=Matriz_estados[f][c];

            switch(matriz_e){

        case 0: {if(Matriz_probabilidades[f][c]<(valor)){Matriz_pasos[f][c]=1;}
                    else {Matriz_pasos[f][c]=0;}

                }break;

        case 1: {if((Matriz_estados[f][sum[c]]==2 || Matriz_estados[f][res[c]]==2 || Matriz_estados[res[f]][c]==2 || Matriz_estados[sum[f]][c]==2) && (Matriz_inmunity[f][c] < (1-(g))))
                   {Matriz_pasos[f][c]=2;}

                /* if((i==1) && (Matriz_estados[sum[f]][c]!=2) && (Matriz_estados[res[f]][c]!=2) && (Matriz_estados[f][sum[c]]!=2) && (Matriz_estados[f][res[c]]!=2) && (Matriz_f[f][c]<h))
                   {Matriz_pasos[f][c]=2;}*/

                }break;

        case 2: {Matriz_pasos[f][c]=0;
                }break;

                                }
                            }
                        }

        for(int f=0;f<fil;f++){
            for(int c=0;c<col;c++){
                    Matriz_estados[f][c] = Matriz_pasos[f][c];
                            }
                        }
                    }

///cierra el archivo
//matriz_moleculas.close();
f_vs_t.close();

}
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/

int main(){

/// generar el patr�n de quemado en la matriz
    srand(time(0)); contorno(); randomizar();Condiciones_Periodicas();

largest_label=0;

Matriz_Moleculas(); Matriz_Estados(); Matriz_Pasos(); Matriz_Occupied();
            //cout<<"L = "<<fil<<endl;
            //cout<<"Pasos = "<<pasos<<endl;
            //cout<<"repeticiones = "<<repeticiones<<endl;
            //cout<<"P = "<<P<<endl<<endl;
            //cout<<"Tama�o    Cantidad"<<endl;

        for(int j=0;j<repeticiones;j++){
           // Matriz_Moleculas(); Matriz_Estados(); Matriz_Pasos(); Matriz_Occupied();

            fuegos=0;
            vacios=0;
            arboles=0;
            fvt=0;
            vvt=0;
            avt=0;

            Pasos(P);
//cout<<endl;
            }
            
      //  }

    return 0;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/

///resetea el histograma DF
/*
    for(int i=0;i<NHS+1;i++){Histograma[i]=0;}
        for (int i = 0; i < NCAJAS; i++) {
            ncajas[i] = 0;
    }
*/
/*
    for(int k = 0; k < NCAJAS; k++) {
        L = LMIN + k * (LMAX - LMIN) / (NCAJAS - 1); // tama�o de la caja
            cout<<L<<"       "<<(ncajas[k])/(pasos-salteo)<<endl;
    }
*/

///imprimir cosas en  pantalla
//for(int i=0;i<NHS+1;i++){cout<<(float)(i*fil)/NHS<<" "<<Histograma[i]<<endl;}cout<<endl<<endl;

/*
 ofstream histograma0("0histograma.txt");// to_string() need c++11
 for(int i=0;i<(NHS+1);i++){histograma0<<(int)(i*fil)/NHS<<" "<<Histograma[i]<<endl;}
 histograma0.close();
*/
/*
 ofstream fuego("densidades.txt");
 densidades<<"fuegos= "<<fuegos<<endl<<"arboles= "<<arboles<<endl<<"vacios= "<<vacios<<endl;
 densidades.close();
                if(Matriz_estados[f][col]==2 || Matriz_estados[f][0]==2 || Matriz_estados[fil][c]==2 || Matriz_estados[0][c]==2)
                    {cout<<i<<endl;
                    Distribucion();}
                    break;
*/
/*
     string base("histograma.txt");
     ofstream histograma(to_string(i)+base);// to_string() need c++11
     for(int i=0;i<(NHS+1);i++){histograma<<(int)(i*fil)/NHS<<" "<<Histograma[i]<<endl;}
     histograma.close();
    }
*/
