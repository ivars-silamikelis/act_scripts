#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#define READ_LENGTH 200
#define REPLICATION_NUMBER 1
void make_read(char * s);
char * char_repeat( int n, char c);
char * create_revcom(char * sequence);
/* no DNS sekvences izveido 200 bp garus nolasījumus ar noteiktu pārklājumu */
//Compile it with gcc shredder4.c -o shredder
int main(int argc, char *argv[])
{ 
  
  char  *sequence,*temp;
  sequence=malloc(sizeof(char));
  int is_sequence=0;
  FILE *filepointer;
  int character;
  int nnuc=0;
  if (argv[1] == NULL){
    printf("Usage:\n");
    printf("./shredder input_filename.fa >output_filename.fastq\n");
    return 1;
  }
  filepointer=fopen(argv[1],"read");
  while ((character=fgetc(filepointer))!=EOF){
    if (is_sequence==1  && character!='\n' && character!=EOF){
      sequence[nnuc]=character;				//sequence ir datu masīvs, kurā ir no fasta faila nolasītā DNS sekvence
      temp=realloc(sequence, (nnuc+2)*sizeof(char));	//temp ir pagaidu datu masīvs, kas ir sequence kopija, tikai tam tiek piešķirta lielāka atmiņa
      if (temp !=NULL){
        sequence=temp;
      }else{
        free(sequence);
        printf("Error allocating memory!\n");
        return 1;
      }
      ++nnuc;
    }
    if (character=='>'){
      is_sequence=0;
    }
    if (character=='\n' && is_sequence==0){
      is_sequence=1;
    }
  }
  sequence[nnuc]='\0';
  
  //Izveidotais datu masīvs ar DNS sekvenci tiek norādīts kā arguments nolasījumu veidojošajai funkcijai make_read
  make_read(sequence);
  
  //Beigās tiek atbrīvota datu masīvam piešķirtā atmiņa 
  free(sequence);
  fclose(filepointer); 
  return 0;
}
/* make_read ir funkcijas, kas nolasa DNS sekvenci un izveido nolasījumus ar noteiktu garumu (200) un intervālu (2) */
void make_read(char * s)
{ 
  
  char * revcom_s=create_revcom(s);				//apgriezti komplementāra dotā DNS sekvence, no kuras tiks veidoti apgriezti komplementāri nolasījumi
  int replication_of_read=REPLICATION_NUMBER;	//Cik reizes viens nolasījums tiks atkārtots
  int n_replicated;  
  int length=0;
  int i;										//pozīcija DNS sekvencē, sākot no kuras tiek veidoti nolasījumi pēc dotās sekvences
  int j;
  int is_revcom=0;								//Norāda veidojamā nolasījuma virzienu
  int id=0;
  int read_size=READ_LENGTH;					//nolasījumu garums
  char read;
  char * quality = char_repeat(read_size, 'M');
  while(s[length]!='\0'){
    length++;
  }
  //Izveido gan tiešus, gan apgriezti komplementārus nolasījumus, katram nukleotīdam piešķirot kvalitāti 34 (ASCII simbols M)
  //Nolasījumi tiek veidoti ik pēc 2 nukleotīdiem. Rezultātā pārklājums būs: nolasījuma garums/nolasījumu intervāls, 200/2=100
  //200 nukleotīdu garie nolasījumi tiek veidoti sākot no pozīcijas 0 DNS sekvencē. Turpmākie nolasījumi tiek veidoti par 2 nukleotīdiem tālāk
  for (i=0;(i+read_size)<length;i+=2){
    for (n_replicated=1;n_replicated<=replication_of_read;++n_replicated){
      printf("@read_%d\n",id);
      ++id;
	  if (is_revcom==0){
        for (j=i;j<=i+read_size-1;++j){
          printf("%c",s[j]);
        }
		is_revcom=1;
      } else if (is_revcom==1){
        for (j=i;j<=i+read_size-1;++j){
          printf("%c",revcom_s[j]);
        }
        is_revcom=0;
    }
    printf("\n+\n");
    puts(quality);

    }
  } 
  free(quality);
  free(revcom_s);
}
/* funkcija kas ļauj atkārtot noteiktu simbolu n reizes */ 
char * char_repeat( int n, char c){
  char * dest = malloc(n+1);
  memset(dest, c, n);
  dest[n]='\0';
  return dest;
}

// funkcija, kas izveido sekvencei apgrieztu komplementāru sekvenci
char * create_revcom(char * sequence )
{
  int i;
  int j=0;
  char * revcom;
  int last_index=strlen(sequence)-1;
  revcom = malloc((last_index+2)*sizeof(char));
  for (i=last_index;i>=0;--i){
    if (sequence[i]=='A'){
      revcom[j]='T';
    } else if  (sequence[i]=='T'){
      revcom[j]='A';
    } else if (sequence[i]=='G'){
      revcom[j]='C';
    } else if (sequence[i]=='C'){
      revcom[j]='G';
    } else if (sequence[i]=='N'){
      revcom[j]='N';
    } else {
      revcom[j]='N';
      fprintf(stderr, "%c nestandarta nukleotīds!\n",revcom[j]);
    }
    ++j;
  }
  revcom[j]='\0';
  return revcom;
}


