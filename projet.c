// gcc  -std=c99 projet.c -lm -O3 -lsndfile -o projet.o

/* carte :   
    source /opt/poky/1.7.1/environment-setup-cortexa5hf-vfp-poky-linux-gnueabi
  $CC -std=c99 projet.c -lm -O3 -lsndfile -o projet.o
  sudo minicom -D /dev/ttyACM0 -b 115200
  mount /dev/sda1 /media
  cd /media
  ./projet.o
*/


#define _XOPEN_SOURCE 500  //pour faire fonctionner msleep

#include <stdio.h>   // pour printf
#include <stdlib.h>  // pour exit 
#include <fcntl.h>   // pour open, read, write et close
#include <unistd.h>  // pour sleep
#include <errno.h>   // pour perror
#include <time.h>
#include <sndfile.h>
#include <complex.h>
#include <math.h>
#include <time.h>



#define windowsize  128 //taille fenetre de spectre
#define SWAP(a,b) ctmp=(a); (a)=(b); (b)=ctmp
#define ARRAY_LEN(x) ((int) (sizeof (x) / sizeof (x [0])))
#define MAX(x,y) ((x) > (y) ? (x) : (y))
#define MIN(x,y) ((x) < (y) ? (x) : (y))
//#define M_PI 3.14159265359


const int datalen=1024; //taille des échantillons testés;
const int H = 30;   //hauteur fenetre

/*-------------------------------Fonctions de sons-----------------------------*/
sf_count_t sfx_mix_mono_read_double(SNDFILE * file, double * data, sf_count_t datalen){
  SF_INFO info;
  static double multi_data [2048] ;
  int k, ch, frames_read ;
  sf_count_t dataout = 0 ;
  sf_command (file, SFC_GET_CURRENT_SF_INFO, &info, sizeof (info)) ;
  if (info.channels == 1)
    return sf_read_double (file, data, datalen) ;
  while (dataout < datalen)
    {   int this_read ;
      this_read = MIN (ARRAY_LEN (multi_data) / info.channels, datalen) ;
      frames_read = sf_readf_double (file, multi_data, this_read) ;
      if (frames_read == 0)
        break ;
      for (k = 0 ; k < frames_read ; k++)
        {       double mix = 0.0 ;
          for (ch = 0 ; ch < info.channels ; ch++)
            mix += multi_data [k * info.channels + ch] ;
          data [dataout + k] = mix / info.channels ;
        } ;
      dataout += this_read ;
    } ;
  return dataout ;
}


/*---------------------------------Fonctions de FFT-----------------------------*/


//Transforme un tableau de double en tableau de complex
complex * array_double_to_complex(double *data, int size, complex *data_complex){
  //On parcourt le tableau de double
  for (int i = 0; i < size; i++){
    //on transforme les double en complex
    data_complex[i] = data[i] + I*0.0; 
  }
  return data_complex;
}

//calcule log2(n)
double log2(double n){
    return log10(n) / log10(2);
}


// Remplis le tableau *result de valeur de FFT de l'entree *data
void fftrec(complex *data, complex *result, unsigned int size, int log2n)
{
  complex ypair[size], yimpair[size], Fimpair[size], Fpair[size];
  int n, k, N2;

  if (size > 1)
    {
    N2 = size/2;
    for ( n=0;n<N2 ; n++)
      {
      ypair[n] = data[n] + data[n+N2];
      yimpair[n] = (data[n] - data[n+N2]) * cexp (-2*I*M_PI*n/size);
    }
    fftrec(ypair, Fpair, N2, log2n);
    fftrec(yimpair, Fimpair, N2, log2n);
    for ( n=0;n<N2 ; n++)
      {
      result[2*n] = Fpair[n];
      result[2*n+1] = Fimpair[n];
    }
  }
  else
    {
    result[0] = data[0];
    return ;
  }
}


//Transforme un tableau de complexe en tableau de module
double *array_complex_to_module(complex *data_complex, int size, double *data_module){
    for (int i = 0; i < size; i++){
        data_module[i] = sqrt(creal(data_complex[i])*creal(data_complex[i]) + cimag(data_complex[i])*cimag(data_complex[i]));
    }
    return data_module;
}


// Cree le spectre de fréquences des modules de la FFT
double *spectre_a_l_echelle(double *data, int size, double *spectrum){
    double moyenne = 0.0;
    int position = 0;
    for (int i = 0; i < size; i++){
        // ici on fait une moyenne sur G valeurs
        moyenne += data[i];
        if (i % (size / windowsize) == 0 && i > 0)
        {
            spectrum[position] = moyenne;
            position += 1;
            moyenne = 0;
        }
    }
    return spectrum;
}


//delai
int msleep(unsigned int tms){
    return usleep(tms * 1000);
}


// trace le spectre mis à l'échelle sur le terminal
void affichage_spectre(double *spectrum, int samplerate){ 
  //on cherche le max de spectrum
  double *p = spectrum;
  double max=*p;
  for(int i=0 ; i<windowsize ; i++){
    max=MAX(max,*(p+1));
    p++;
  }

  int new_spectrum[windowsize];   //le spectre mis à l'échelle de la fenetre
  for (int i = 0; i < windowsize; i++)
        /*on met le spectre à l'échelle de la fenetre en fonction du max trouvé*/

        /* max evolutif*/
        //new_spectrum[i] = spectrum[i] * H / (max+max/10);

        /* max fixé */
        new_spectrum[i] = spectrum[i] * H / (80+80/10);


  //affichage du spectre (en couleur) de haut en bas
  for (int j = H; j > 0 ; j--){
    if( j>3*H/4 ){
      printf("\033[32m"); //couleur
      for (int i = 0; i < windowsize; i++)
          (j > new_spectrum[i]) ? printf(" ") : printf("8");
      printf("\n");
    }

     if( (j>H/2) && (j<3*H/4) ){
      printf("\033[33m"); //couleur
      for (int i = 0; i < windowsize; i++)
          (j > new_spectrum[i]) ? printf(" ") : printf("8");
      printf("\n");
    }

     if( (j>H/4) && (j<H/2) ){
      printf("\033[31m"); //couleur
      for (int i = 0; i < windowsize; i++)
          (j > new_spectrum[i]) ? printf(" ") : printf("8");
      printf("\n");
    }

     if( (j<H/4) ){
      printf("\033[35m"); //couleur
      for (int i = 0; i < windowsize; i++)
          (j > new_spectrum[i]) ? printf(" ") : printf("8");
      printf("\n");
    }
  }
    
  //barre de limite horizontale avant le reflet
  printf("\033[34m"); //couleur
  for (int i = 0; i < windowsize; i++){
    printf("8");
  }
  printf("\n");

  //reflet
  for (int j = 1; j < H/4 ; j++){
    if( j<H/4 ){
      printf("\033[35m"); //couleur
      for (int i = 0; i < windowsize; i++)
          (j > new_spectrum[i]) ? printf(" ") : printf("8");
      printf("\n");
    }
  }



  printf("\033[37m");  //couleur
  // barre d'echelle des frequences
  for (int i = 0; i < windowsize; i++)
      printf("_");
  printf("\n0 Hz");
  for (int i = 0; i < windowsize - 11; i++)
      printf(" ");
  printf("%d Hz\n", samplerate / 4);
}

/*------------------------------------main----------------------------------------*/
int main(void) {
    SNDFILE *infile;
    SF_INFO sfinfo;

    //fichier son
    char input[]="Sons/Video/dance.wav";


    //ouverture du fichier son
    infile = sf_open (input, SFM_READ, &sfinfo);
    if(infile == NULL){
      /* Open failed so print an error message. */
      printf("Not able to open input file %s.\n",input);
      /* Print the error message fron libsndfile. */
      sf_perror(NULL) ;
      return 1 ;
    }

    //definition des variables d'analyse de son
    int samplerate = sfinfo.samplerate;
    int channels = sfinfo.channels;
    //int datalen=1024; //taille des échantillons testés;
    double data[datalen];
    int readcount;
    // definition des variables de fft
    complex data_complex[datalen];
    complex data_complex_fft[datalen];
    double data_module[datalen];
    double spectrum[windowsize];
    clock_t previous, diff;
    double msec;
    int delai = datalen * 1000 / samplerate;   //temps d'attente

    /*boucle d'analyse de frequences*/
    while ((readcount = sfx_mix_mono_read_double (infile, data, datalen)) > 0){
        //FFT à faire dans ce while//
        previous = clock();
        //on transforme la data du son en complexe
        array_double_to_complex(data, datalen, data_complex);
        //on calcule le FFT sur la data complexe
        fftrec(data_complex, data_complex_fft, datalen, log2(datalen));
        //on transforme la data complexe de la fft en double
        array_complex_to_module(data_complex_fft, datalen, data_module);
        
        //cree les valeurs du spectre
        spectre_a_l_echelle(data_module, datalen/4, spectrum);
        //trace le spectre sur le terminal 
        affichage_spectre(spectrum, samplerate);
        
        //gere le temps d'affichage du spectre a l'ecran, et fais un clear
        diff = clock() - previous;
        msec = delai - (diff * 1000 / CLOCKS_PER_SEC);
        msleep(msec);
        //msleep(80);
        printf("\e[1;1H\e[2J");
    

    }

    /* ferme le fichier son */
    sf_close(infile) ;
  return EXIT_SUCCESS;
}
