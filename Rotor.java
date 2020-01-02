public class Rotor{
   int[] rotor;
   int turn;
   int signal;
   
   //Each index holds the position of a number
   public Rotor(int start_rot, int set){
      this.turn = start_rot;
      //this.rotor = new int[]{18,17,16,12,21,6,19,13,14,0,20,5,4,11,2,10,1,9,8,24,25,22,3,23,7,15};
      if(set == 1){
         this.rotor = new int[]{18,17,16,12,21,6,19,13,14,0,20,5,4,11,2,10,1,9,8,24,25,22,3,23,7,15};
         this.signal = 24;
      }
      else if(set == 2){
         this.rotor = new int[]{11,15,16,5,1,14,6,9,24,19,18,10,20,21,17,2,0,25,7,8,3,23,4,22,12,13};
         this.signal = 5;
      }
      else if(set == 3){
         this.rotor = new int[]{25,23,19,2,6,4,18,8,1,5,16,7,15,14,22,11,17,12,0,10,21,3,9,20,24,13};
         this.signal = 14;
      }
      else if(set == 4){
         this.rotor = new int[]{8,0,11,1,7,25,6,16,24,12,22,9,3,21,13,18,17,10,5,4,2,15,20,14,19,23};
         this.signal = 13;
      }
      else if(set == 5){
         this.rotor = new int[]{4,8,0,19,21,1,3,7,14,6,9,15,16,22,11,5,24,20,23,18,17,2,25,12,13,10};
         this.signal = 9;
      }
      
   }
   
   public int forwardTranslate(int n){
      return this.rotor[(n+this.turn)%26];
   }
   
   public int backTranslate(int n){
      int i = 0;
      int out = 0;
      while(this.rotor[i]!=n){
         i++;
      }
      out = (i-this.turn)%26;
      if(out<0) return out+26;
      else return out;
   }
   
   public void rotate(){
      turn++;
   }
   
}