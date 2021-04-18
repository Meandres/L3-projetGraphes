#include <cstdlib>
#include <iostream>
#include <cmath>
#include <new>
#include <fstream>
#include <time.h>
#include <vector>
#include <random>
#include <cstring>

using namespace std;

long n; // nombre de sommets
long N; //nb n des graphes de petersen généralisés
long **adj; //[n][n];  // matrice d'adjacence du graphe
long *couleur1; //[n];  // couleurs des sommets pour l'agorithme exact
long *couleur2; //[n]; // couleurs pour DSATUR
long *couleurTamp;
long *DSAT; //[n]; // degr�s de saturation
long *Degre; //[n]; // degr�s des sommets
bool trouve=false; // permet de stopper l'algorithme exact
                   // quand une coloration  a ete trouvee

int H=2, K=1; //paramètres de la coloration

void genereGP(long k){
  for(long i=0; i<N; ++i){
    for(long j=0; j<N; j++){
      adj[i][j]=0; adj[i+N][j]=0; adj[i][j+N]=0; adj[i+N][j+N]=0;
      if(j==(i+1)%N || j==(i-1)%N){
        adj[i][j]=1;
      }
      if(i==j){
        adj[i+N][j]=1;
        adj[i][j+N]=1;
      }
      if(j==(i+k)%N || j==(i-k)%N){
        adj[i+N][j+N]=1;
      }
    }
  }
}
void genereG(){
  srand(time(NULL));
  vector<long> pos;
  long rd, i;
  for(i=0; i<N; ++i){
    for(long j=0; j<N; ++j){
      adj[i][j]=0; adj[i+N][j]=0; adj[i][j+N]=0; adj[i+N][j+N]=0;
    }
    pos.push_back(N+i);
  }
  for(i=0; i<N; ++i){
    adj[i][(i+1)%N]=1;
    if(i==0)
      adj[i][N-1]=1;
    else
      adj[i][(i-1)%N]=1;

    adj[i+N][(i+1)%N + N]=1;
    if(i==0)
      adj[i+N][N-1 + N]=1;
    else
      adj[i+N][(i-1)%N + N]=1;

    rd=rand()%pos.size();
    adj[i][pos.at(rd)]=1;
    adj[pos.at(rd)][i]=1;
    pos.erase(pos.begin()+rd);
  }
}
bool convientL21(long x, long c) // teste si la couleur c peut �tre donnee au sommet x
{
     for(long i=0;i<x;i++)
      if(adj[x][i]){
        if(abs(couleur1[i] - c)<H) //test si la couleur des voisins est éloignée d'au moins H
          return false;
        for(long j=0; j<n; ++j)
          if(adj[i][j] && abs(couleur1[j] - c)<K) return false; //test si la couleur des voisins des voisins est éloignée d'au moins K
      }

     return true;
}
bool convientDSATL21(long x, long c) // teste si la couleur c peut �tre donnee au sommet x
{
     for(long i=0;i<n;i++)
      if(adj[x][i]){
        if(abs(couleur2[i] - c)<H) //test si la couleur des voisins est éloignée d'au moins H
          return false;
        for(long j=0; j<n; ++j)
          if(adj[i][j] && abs(couleur2[j] - c)<K) return false; //test si la couleur des voisins des voisins est éloignée d'au moins K
      }

     return true;
}

void colorRR(long x, long k) // fonction recursive pour tester toutes les couleurs possible pour le sommet x
{
     if(x==n)
     { //cout << "Coloration en " << k << " couleurs trouvee" << endl;
       //for(int i=0;i<n;i++) cout << "couleur de " << i << " : " << couleur1[i] << endl; //int z;cin >> z;
       trouve=true;
       for (long i = 0; i < N*2; i++) {
         couleurTamp[i] = couleur1[i];
       }
     }
     else
     for(long c=1;c<=k;c++)
      if(convientL21(x,c))
	{couleur1[x]=c;//cout << "=>couleur de " << x << " : " << couleur[x] << endl;
         colorRR(x+1,k);
	 if(trouve) return;}
     //return false;
}


void colorexact(long k) // teste si le graphe possede une coloration en k couleurs en essayant toutes les combinaisons
{
    for(long i=0;i<n;i++){
      couleur1[i]=0;
    }
    colorRR(0,k);
     //if(!trouve) cout << "Pas de coloration en " << k <<" couleurs" << endl;
}



long nbChromatique(long d) // calcule le nombre chromatique en testant à partir de d couleurs et diminuant k tant que c'est possible
{
  long k=d+1;
  do {
      k--;
      trouve=false;
      colorexact(k);
     }
  while(trouve);
  return k+1;
}


long dsatMax()
{
  long maxDeg=-1,maxDSAT=-1,smax=0;
  for(long i=0;i<n;i++)
  if(couleur2[i]==0 && (DSAT[i]>maxDSAT || (DSAT[i]==maxDSAT && Degre[i]>maxDeg)))
   { maxDSAT=DSAT[i]; maxDeg=Degre[i]; smax=i;}
  return smax;
}

long DSATUR(){
  long nb=0,c,x,cmax=0;
  for(long i=0;i<n;i++)
  {
    couleur2[i]=0; DSAT[i]=0; Degre[i]=0;
    for(long j=0;j<n;j++)
     if(adj[i][j]) Degre[i]++;
    DSAT[i]=Degre[i];
  }

  while(nb<n)  // tant qu'on a pas colori� tous les sommets
  {
    c=1;
    x=dsatMax(); // on choisit le sommet de DSAT max non encore colori�
    while(!convientDSATL21(x,c)) c++; // on cherche la plus petite couleur disponible pour ce sommet
    for(long j=0; j<n;j++) // mise � jour des DSAT
     {
       if(adj[x][j] && convientDSATL21(j,c)) DSAT[j]++; // j n'avait aucun voisin colori� c,on incr�mente donc son DSAT
     }
    couleur2[x]=c;
    if(cmax<c) cmax=c;
    nb++;
  }
  return cmax;
}


void affichageDSAT(long k, long nbc){
    cout << "n : " << N << " k : " << k << endl;
    cout << "Nb chromatique :" << nbc << endl;
    for (long i = 0; i < N ; i++) {
        cout << couleur2[i] << " ";
    }
    cout << endl;
    for (long i = N; i < N*2 ; i++) {
        cout << couleur2[i] << " ";
    }
    cout << endl << endl;
}

void affichageColorExact(long k, long nbc){
    cout << "n : " << N << " k : " << k << endl;
    cout << "Nb chromatique :" << nbc << endl;
    for (long i = 0; i < N ; i++) {
        cout << couleurTamp[i] << " ";
    }
    cout << endl;
    for (long i = N; i < N*2 ; i++) {
        cout << couleurTamp[i] << " ";
    }
    cout << endl << endl;
}

void outputTempsExec(string fic, bool prismes){
  float tempsColorExact, tempsDSAT;
  long k;
  clock_t t1, t2;
  ofstream myfile;
  myfile.open(fic);
  myfile << "Valeur de n; Valeur de k; temps exec colorExact\n";
  for (int i = 10 ; i <= 20 ; i+=2) {
    for (int j = 2; j <= 4; j++) {
      cout << i << " " << j << endl;
      N = i;
      k=floor(N*j/6);
      n=2*N;
      adj=new long*[n];
      for (long x = 0; x < n; x++)
        adj[x] = new long[n];
      couleur1= new long[n]; couleurTamp = new long[n];
      if(prismes)
        genereG();
      else
        genereGP(k);

      t1 = clock();
      nbChromatique(N);
      t2 = clock();
      tempsColorExact = (double)(t2 - t1)/CLOCKS_PER_SEC;

      free(adj); free(couleur1); free(couleurTamp);

      myfile << N << ";" << k << ";" << tempsColorExact << "\n";
    }
  }
  cout << "dsat" << endl;
  myfile << "Valeur de n; Valeur de k; temps exec DSAT\n";
  for(int i=8; i<=13; ++i){
    for(int j=2; j<=4; ++j){
      cout << i << " " << j << endl;
      N = pow(2, i);
      k=floor(N*j/6);
      n=2*N;
      adj=new long*[n];
      for (long x = 0; x < n; x++)
        adj[x] = new long[n];
      couleur2 = new long[n]; couleurTamp = new long[n];
      DSAT = new long[n]; Degre = new long[n];

      if(prismes)
        genereG();
      else
        genereGP(k);

      t1 = clock();
      DSATUR();
      t2 = clock();
      tempsDSAT = (double)(t2 - t1)/CLOCKS_PER_SEC;

      free(couleur2); free(couleurTamp); free(DSAT); free(Degre);

      myfile << N << ";" << k << ";" << tempsDSAT << "\n";
    }
  }
  myfile.close();
}


int main(int argc, char *argv[])
{
  bool correct=true, optG=false, optO=false;
  long p,k,nbc;
  string opt;
  if(argc<=2){
    correct=false;
  }
  if(correct){
    if(strcmp(argv[1], "-o")==0){
      if(argc>=3){
        optO=true;
        opt=argv[2];
      }
      else{
        correct=false;
      }
      if(argc==4){
        if(strcmp(argv[3], "-g")==0){
          optG=true;
        }
      }
    }
    else{
      N = atoi(argv[1]);
      if(argc==3){
        if(strcmp(argv[2], "-g")==0)
          optG=true;
        else
          k=atoi(argv[2]);
      }
      else
        correct=false;
    }
  }


if(!correct){
  cerr << "Usage : main N K/-g\nOu : main -o nomFic (-g)\nAvec N le nombre de sommet d'un cycle.\nK le parametre des graphes de Petersen ou -g qui indique que l'on veut les prismes.\n-o nomFic qui permet la sortie dans un fichier csv" << endl;
  exit(1);
}

  if(optO){
    outputTempsExec(opt, optG);
  }
  else{
    n=2*N;
    adj=new long*[n];
    for (long i = 0; i < n; i++)
       adj[i] = new long[n];
    couleur1= new long[n]; couleur2 = new long[n]; couleurTamp = new long[n];
    DSAT = new long[n]; Degre = new long[n];

    if(optG)
      genereG();
    else
      genereGP(k);

    cout << " DSATUR " << endl;
    nbc=DSATUR();
    affichageDSAT(k, nbc);

    cout << "ColorExact :" << endl;
    nbc=nbChromatique(n);
    affichageColorExact(k, nbc);
  }

  return 0;

}
