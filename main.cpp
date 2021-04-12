#include <cstdlib>
#include <iostream>
#include <cmath>
#include <new>

using namespace std;

int n; // nombre de sommets
int N; //nb n des graphes de petersen généralisés
int **adj; //[n][n];  // matrice d'adjacence du graphe
int *couleur1; //[n];  // couleurs des sommets pour l'agorithme exact
int *couleur2; //[n]; // couleurs pour DSATUR
int *couleurTamp;
int *DSAT; //[n]; // degr�s de saturation
int *Degre; //[n]; // degr�s des sommets
bool trouve=false; // permet de stopper l'algorithme exact
                   // quand une coloration  a ete trouvee

int H=2, K=1; //paramètres de la coloration

void genereGP(int k){
  for(int i=0; i<N; ++i){
    for(int j=0; j<N; j++){
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
  vector<int> pos;
  int rd, i;
  for(i=0; i<N; ++i){
    for(int j=0; j<N; ++j){
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
bool convientL21(int x, int c) // teste si la couleur c peut �tre donnee au sommet x
{
     for(int i=0;i<x;i++)
      if(adj[x][i]){
        if(abs(couleur1[i] - c)<H) //test si la couleur des voisins est éloignée d'au moins H
          return false;
        for(int j=0; j<n; ++j)
          if(adj[i][j] && abs(couleur1[j] - c)<K) return false; //test si la couleur des voisins des voisins est éloignée d'au moins K
      }

     return true;
}
bool convientDSATL21(int x, int c) // teste si la couleur c peut �tre donnee au sommet x
{
     for(int i=0;i<n;i++)
      if(adj[x][i]){
        if(abs(couleur2[i] - c)<H) //test si la couleur des voisins est éloignée d'au moins H
          return false;
        for(int j=0; j<n; ++j)
          if(adj[i][j] && abs(couleur2[j] - c)<K) return false; //test si la couleur des voisins des voisins est éloignée d'au moins K
      }

     return true;
}

void colorRR(int x, int k) // fonction recursive pour tester toutes les couleurs possible pour le sommet x
{
     if(x==n)
     { //cout << "Coloration en " << k << " couleurs trouvee" << endl;
       //for(int i=0;i<n;i++) cout << "couleur de " << i << " : " << couleur1[i] << endl; //int z;cin >> z;
       trouve=true;
       for (int i = 0; i < N*2; i++) {
         couleurTamp[i] = couleur1[i];
       }
     }
     else
     for(int c=1;c<=k;c++)
      if(convientL21(x,c))
	{couleur1[x]=c;//cout << "=>couleur de " << x << " : " << couleur[x] << endl;
         colorRR(x+1,k);
	 if(trouve) return;}
     //return false;
}


void colorexact(int k) // teste si le graphe possede une coloration en k couleurs en essayant toutes les combinaisons
{
    for(int i=0;i<n;i++)
     couleur1[i]=0;
     colorRR(0,k);
     //if(!trouve) cout << "Pas de coloration en " << k <<" couleurs" << endl;
}



int nbChromatique(int d) // calcule le nombre chromatique en testant � partir de d couleurs et diminuant k tant que c'est possible
{
  int k=d+1;
  do {
      k--;
      trouve=false;
      colorexact(k);
     }
  while(trouve);
  return k+1;
}


int dsatMax()
{
  int maxDeg=-1,maxDSAT=-1,smax=0;
  for(int i=0;i<n;i++)
  if(couleur2[i]==0 && (DSAT[i]>maxDSAT || (DSAT[i]==maxDSAT && Degre[i]>maxDeg)))
   { maxDSAT=DSAT[i]; maxDeg=Degre[i]; smax=i;}
  return smax;
}

int DSATUR()
{
  int nb=0,c,x,cmax=0;
  for(int i=0;i<n;i++)
  {
    couleur2[i]=0; DSAT[i]=0; Degre[i]=0;
    for(int j=0;j<n;j++)
     if(adj[i][j]) Degre[i]++;
    DSAT[i]=Degre[i];
  }

  while(nb<n)  // tant qu'on a pas colori� tous les sommets
  {
    c=1;
    x=dsatMax(); // on choisit le sommet de DSAT max non encore colori�
    while(!convientDSATL21(x,c)) c++; // on cherche la plus petite couleur disponible pour ce sommet
    for(int j=0; j<n;j++) // mise � jour des DSAT
     {
       if(adj[x][j] && convientDSATL21(j,c)) DSAT[j]++; // j n'avait aucun voisin colori� c,on incr�mente donc son DSAT
     }
    couleur2[x]=c;
    if(cmax<c) cmax=c;
    nb++;
  }

  return cmax;
}
void affichageDSAT(int k){
    cout << "n : " << N << " k : " << k << endl;
    for (int i = 0; i < N ; i++) {
        cout << couleur2[i] << " ";
    }
    cout << endl;
    for (int i = N; i < N*2 ; i++) {
        cout << couleur2[i] << " ";
    }
    cout << endl <<endl;
}

void affichageColorExact(int k){
    cout << "n : " << N << " k : " << k << endl;
    for (int i = 0; i < N ; i++) {
        cout << couleurTamp[i] << " ";
    }
    cout << endl;
    for (int i = N; i < N*2 ; i++) {
        cout << couleurTamp[i] << " ";
    }
    cout << endl <<endl;
}


int main()
{
  int p,k,K,nbc;
  cout << "nombre N des graphes de petersen généralisés (=> nb sommet/2)" << endl;
  cin >> N;
  //cout << "proba d'arete: " << endl;
  //cin >> p;
  cout << "nombre k" << endl;
  cin >> K;

  n=2*N;
  adj=new int*[n];
  for (int i = 0; i < n; i++)
     adj[i] = new int[n];
  couleur1= new int[n]; couleur2 = new int[n]; couleurTamp = new int[n];
  DSAT = new int[n]; Degre = new int[n];

  genereGP(k);
  //genereG();

  cout << " DSATUR " << endl;
  k=DSATUR();
  affichageDSAT(K);



  cout << "ColorExact :" << endl;
  nbc=nbChromatique(n);
  affichageColorExact(K);

  return 0;
}
