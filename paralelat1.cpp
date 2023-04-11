#include <iostream>
#include <cmath>
#include <cctype>
#include <vector>
#include <string>
#include <algorithm>
using namespace std;
/*
void ferrari(int gradomaximo, vector cofeciente2[5]) // intento de calcular las raices para una ecuacion cuartica
{
    if(gradomaximo != 4)
    {
        cout<< "No corresponde a este metodo" <<endl;
    }
    double aux1, aux2, p,q,r, pp, rr, qq, x, t, t2, t3, t4; //los aux son para calcular las raices y p,q y r son parte de la formula
    p= ((8*cofeciente2[1])-(3*pow(cofeciente2[0],2)))/8;
    q= ((8*cofeciente2[2])-(4*cofeciente2[0]*cofeciente2[1])+(pow(cofeciente2[0],3)))/8;
    r= ((256*cofeciente2[3])-(64*cofeciente2[0]*cofeciente2[2])+(16*cofeciente2[0]*cofeciente2[0]*cofeciente2[1])-(3*pow(cofeciente2[0],4)))/256;
    pp= cbrt((-p/4)+((pow(p/4,2))+ sqrt(pow((pow(p,3)/256),3)))) + cbrt((-p/4)- (pow(p/4,2))+ (pow((pow(p,3)/256),3))); //P = [(-p/4) + (((p/4)^2) + (p^3/256))^3)^(1/2))^(1/3)] + [(-p/4) - (((p/4)^2) + (p^3/256))^3)^(1/2))^(1/3)]
    cout << "La primera raiz es: " << pp << endl;
    qq= sqrt((2*pp-p));
    rr= pp - sqrt(r);
    if(q==0)
    {
        t=(-p+ sqrt(pow(p,2)-r))/2;
        t2=(-p- sqrt(pow(p,2)-r))/2;
        cout << "La primera raiz es: " << t << endl;
        cout << "La segunda raiz es: " << t2 << endl;
    }
    else
    {
    aux1 = (-p + sqrt((p*p)-r))/2;
    aux2 = (-p - sqrt((p*p)-r))/2;
    t= sqrt(aux1)-((cofeciente2[1])/4);
    t2= -(sqrt(aux1))-((cofeciente2[1])/4);
    t3= sqrt(aux2)-((cofeciente2[1])/4);
    t4= -(sqrt(aux2))-((cofeciente2[1])/4);
    cout << "La primera raiz es: " << t << endl;
    cout << "La segunda raiz es: " << t2 << endl;
    cout << "La tercera raiz es: " << t3 << endl;
    cout << "La cuarta raiz es: " << t4 << endl;
    }
} // no resulto por tiempo por ende no se comento que hace 
*/
int gradomaximo(string polinomio){ //obtener el grado mayor y de terminar a que polinio corresponde
 int gradomaximo = 0; //se inicializa el grado mayor en 0
    for (int i = 0; i < polinomio.length(); i++) {
      int exponente = 0;
       if (polinomio[i] == '*') {
            int j = i + 2;
            while (isdigit(polinomio[j])) {//isdigit utilizado para recoorrer el valor del exponente mientras se encuentre un valor numerico a su lado(cuando el exponente es de mas de dos digitos)
                exponente = exponente *10 + (polinomio[j] - '0');//se resta el valor de 0 en ascii al valor obtenido en el polinomio[j], lo que lleva a su correspondiente valor numerico
                j++;
            } 
        } else if (polinomio[i]=='x'&&polinomio[i+1]!='*'){
            exponente =1;
          }
      if (exponente > gradomaximo) {//si se encuentra un grado mayor, se sustituye el valor del grado mayor del polinomio
                gradomaximo = exponente;
            }
    }
return gradomaximo;   //se retorna variable con grado maximo de polinomio
}


vector<string> separarterminos(string polinomio) {//separa cada termino y lo pone en un vector
    vector<string> terminos;//se crea un vector que contendra los terminos por separado del polinomio 
    string termino = "";//variable para ir guardando cada termino del polinomio
    for (int i = 0; i < polinomio.length(); i++) {
        if (polinomio[i] == '+'||polinomio[i] == '-') { // pregunta si es que hya un + o -
            if (!termino.empty()) {
                terminos.push_back(termino);
            }
            termino = "";
            termino += polinomio[i];
        } else if (isdigit(polinomio[i])) { // si es un digito los optines hasta que el caracter sea distinto de un digito
            int j = i;
            while (isdigit(polinomio[j]) ) {
                termino += polinomio[j];
                j++;
            }
            i = j - 1;
        } else if (polinomio[i] == 'x') { // si el caracter es una x 
            termino += polinomio[i];
            if (i < polinomio.length() - 2 && polinomio[i+1] == '*' && polinomio[i+2] == '*') { // comprueba que los siguientes caracteres sean **
                termino += "**";
                i += 2; // se los salta en la proxima iteracion 
            }
        } else if (polinomio[i] == '*') {
            termino += polinomio[i];
        } 
    }
    if (!termino.empty()) {
        terminos.push_back(termino);
    }
    return terminos;
}




vector<double> obtenerCoeficientes(vector<string> terminos) { // obtener los diferentes coeficientes del polinomio
    vector<double> coeficientes;
    for (string termino : terminos) { //recorre la cadena de caracteres
        double coeficiente = 0;
        size_t encuentra = termino.find("x"); // busca en que espacio se encuntra x
        if (encuentra != string::npos) {
            string strCoeficiente = termino.substr(0, encuentra);
            if (strCoeficiente.empty() || strCoeficiente == "+") { // pregunta si esta vacio o un +
                coeficiente = 1;
            } else if (strCoeficiente == "-") { // cpnsulta si hay un -
                coeficiente = -1;
            } else {
                coeficiente = stod(strCoeficiente);
            }
        } else {
            coeficiente = stod(termino);
        }
        coeficientes.push_back(coeficiente);
    }
    return coeficientes;
}








std::vector<int> grados(string polinomio){ // obtienen los grados y los ordena de mayor a menor
  vector<int>grados;
for (int i = 0; i < polinomio.length(); i++) {
        if (polinomio[i]=='x'&& polinomio[i+1] == '*') { //busca cuando este x y el caracter siguiente sea *
            int j = i + 3;
            int exponente = 0;
            while (isdigit(polinomio[j])) { //mientras sea una
                exponente = exponente *10 + (polinomio[j] - '0');
                j++;
            }
          grados.push_back(exponente);
        }
  if (polinomio[i]=='x'&& polinomio[i+1] != '*') {
            int exponente = 1;
          grados.push_back(exponente);
        }
    }
 
   std::vector<int> vectorOrdenado = grados;
    std::sort(vectorOrdenado.begin(), vectorOrdenado.end(), std::greater<int>());
    return vectorOrdenado;
}




std::vector<double> coeficientes(string polinomio,vector<double> coeficiente,int gradomax){// ordena los coeficientes
  
vector<double> grados;
  
for (int i = 0; i < polinomio.length(); i++) {
        if (polinomio[i]=='x'&& polinomio[i+1] == '*') {
            int j = i + 3;
            int exponente = 0;
            while (isdigit(polinomio[j])) { //mientras sea un digito ir guardandoloes en exponente
                exponente = exponente *10 + (polinomio[j] - '0');
                j++;
            }
          
          grados.push_back(exponente);
        }
  if (polinomio[i]=='x'&& polinomio[i+1] != '*') {//cuando este la x sin un numero a su lado se guarde un 1
            int exponente = 1;
          
            grados.push_back(exponente);
          }
        }
   
    int i,j,aux,aux2,min;

    for(i=0;i<grados.size()+1;i++){//usa el ordenamiento por seleccion para ordenar los coeficientes a la vez que se ordenan los grados por orden
      min=i;
      for(j=i+1;j<grados.size()+1;j++){
        if(grados[j]<grados[min]){
            min=j;  
        }
    }
    aux=grados[i];
    aux2=coeficiente[i];
    grados[i]=grados[min];
    coeficiente[i]=coeficiente[min];
    grados[min]=aux;
    coeficiente[min]=aux2;
    }
  
std::reverse(coeficiente.begin(), coeficiente.end());//invierte el orden del vector
 
    return coeficiente;
}




/*
std::vector<double> ordencoef(vector<double> coeficienteordenado, vector<int> gradord,int gradomaximo){ //ordena los coeficientes de acuerdo a su grado


int gradoactual = gradomaximo;
int m=0;
    vector<double> orden(gradomaximo, 0);
     vector<int> totalgrados(gradomaximo, 0);
    for(int j = 0; j < coeficienteordenado.size(); j++){//ordena los coeficientes dependiendo de el vector de grados y el grado máximo

      
        if(gradord[m] == gradoactual){
            totalgrados[j] = gradord[m];
            orden[j]=coeficienteordenado[m];
            gradoactual--;
            m++;
        }else{
            totalgrados[j] = 0;
            gradoactual--;
            orden[j]=0;
        }
    } 

    if(coeficienteordenado.size()>gradord.size()){//agrega el coeficiente que no depende de un grado si la cantidad de coeficientes es mayor a la de grados
      orden.push_back(coeficienteordenado[coeficienteordenado.size()-1]);
     
    }else{
      orden.push_back(0);
    }

  return orden;
}
*/



std::vector<int> obtenergrados(string polinomio){ //obtiene cada grado en el polinomio
  vector<int>grados;
for (int i = 0; i < polinomio.length(); i++) {
        if (polinomio[i]=='x'&& polinomio[i+1] == '*') { //caundo una la incognita x este acompañada de un *
            int j = i + 3;
            int exponente = 0;
            while (isdigit(polinomio[j])) { //buscar mientras sea una digito
                exponente = exponente *10 + (polinomio[j] - '0'); //guardar el exponente
                j++;
            }
          grados.push_back(exponente);
        }
  if (polinomio[i]=='x'&& polinomio[i+1] != '*') { //si no tinene el *, el exponenete es 1
            int exponente = 1;
          grados.push_back(exponente);
        }
    }
 

    return grados;
}



vector<double> raicesPolinomioGrado3(string polinomio) { // calcular raices de polinomio de tercer grado
    // Obtenemos los coeficientes del polinomio
    vector<string> terminos = separarterminos(polinomio);
    vector<double> coeficientes = obtenerCoeficientes(terminos);
    int n = coeficientes.size();
    
    // Obtenemos los grados de las variables
    vector<int> grados = obtenergrados(polinomio);

    // Verificamos que el grado del polinomio sea 3
    if (n != 4 || grados[0] != 3) {
        cout << "El polinomio no es de grado 3 o esta en una forma incompleta" << endl;
        return {};
    }

    //se ingresan los coeficientes
    double a = coeficientes[0];
    double b = coeficientes[1];
    double c = coeficientes[2];
    double d = coeficientes[3];
    double p = (3*a*c - b*b)/(3*a*a);
    double q = (2*b*b*b - 9*a*b*c + 27*a*a*d)/(27*a*a*a);

    // Calculamos las raíces del polinomio utilizando la fórmula de Cardano-Tartaglia
    vector<double> raices;
    double delta = q*q/4 + p*p*p/27;//se ingresa variable delta de la ecuacion
    if (delta > 0) {
        double u = cbrt(-q/2 + sqrt(delta));
        double v = cbrt(-q/2 - sqrt(delta));
        raices.push_back(u + v - b/(3*a));
    } else if (delta == 0) {
        double u = cbrt(-q/2);
        raices.push_back(2*u - b/(3*a));
        raices.push_back(-u - b/(3*a));
    } else {
        double u = 2*sqrt(-p/3);
        double v = acos(-sqrt(-27/(p*p*p))*q/2)/3;
        raices.push_back(u*cos(v) - b/(3*a));
        raices.push_back(u*cos(v + 2*M_PI/3) - b/(3*a));//se utiliza M_PI para obtener un valor mas acertado a pi, se intento con 3.14 pero daba valores ligeramente disntintos y con M_PI es mas preciso
        raices.push_back(u*cos(v + 4*M_PI/3) - b/(3*a));
    }

    return raices;
}




 int main(int argc, char* argv[]) {
    if (argc != 2) {
        cout << "Por favor ingrese el polinomio como argumento al ejecutar el programa" << endl;
        return 0;
    }

    string polinomio = argv[1];
      cout << "Ingrese el polinomio con x como incognita: ";
    getline(cin, polinomio);
  //se ingrasa y se guarda el polinomio
    int gmax = gradomaximo(polinomio);
       cout << "El grado del polinomio es: " << gmax << endl;
    vector<string> terminos = separarterminos(polinomio); // separa los terminos del polinomio
  vector<double>coeficiente=obtenerCoeficientes(terminos); //separa los coeficientes
vector<int> todosgrados =grados(polinomio);
vector<double> coeficientesord =coeficientes( polinomio,coeficiente, gmax);
//vector<double> coeficientesordenados = ordencoef(coeficientesord,todosgrados,gmax);

  
  if(gmax==0){ 
    cout << "Polinomio invalido, reintente "; // en caso que solo se ingresaran numeros se pide reintentarlo
  }
  if(gmax==2){ // en caso de que el polinomio sea una ecuacion cuadratica calcula las raices
    double a, b, c;
    a=coeficiente[0];
    b=coeficiente[1];
    c=coeficiente[2];
double discriminante = b * b - 4 * a * c;
    if (discriminante > 0) {
        double x1 = (-b + sqrt(discriminante)) / (2 * a);
        double x2 = (-b - sqrt(discriminante)) / (2 * a);
        cout << "Las raices son x1 = " << x1 << " y x2 = " << x2 << endl;
    } else if (discriminante == 0) {
        double x = -b / (2 * a);
        cout << "La unica raiz es x = " << x << endl;
    } else {
        double partereal = -b / (2 * a);
        double parteimag = sqrt(-discriminante) / (2 * a);
        cout << "Las raices son x1 = " << partereal << " + " << parteimag << "i y x2 = " << partereal << " - " << parteimag << "i" << endl;
    }

    } 

if(gmax==3){ // en caso de que el polinomio se una ecuacion cubica
       vector<double> raices = raicesPolinomioGrado3(polinomio); // se llama a la funcion para calcular sus raices
    for (double raiz : raices) {
        cout << "Raiz: " << raiz << endl;
    }
  }



   cout <<"integrantes:" << endl;
   cout <<"Mauricio Estrada" << endl;
   cout <<"Oscar Soto" << endl;
   cout <<"Rodrigo Vasquez" << endl;
  return 0;
  }
   

 






