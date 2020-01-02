import java.util.Scanner;

public class MassEnigmaMachine{

   public static void main(String args[]){
        
      Rotor reflector = new Rotor(0,1);
      reflector.rotor = new int[]{25,24,23,22,21,20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0};
      
      Scanner read = new Scanner(System.in);     
      
      Rotor r1 = new Rotor(Integer.valueOf(1), Integer.valueOf(1));
      Rotor r2 = new Rotor(Integer.valueOf(2), Integer.valueOf(2));
      Rotor r3 = new Rotor(Integer.valueOf(3), Integer.valueOf(3));

      
         //message prep
         String test = StdIn.readAll();
         test = test.toLowerCase();
         test = test.replaceAll(" ","");
         test = test.replaceAll("\\p{P}","");
         test = test.replaceAll("\\s+", "");   
              
         //enigma machine
         String output = "";
         for(int i = 0; i<test.length();i++){
            char c = test.charAt(i);
            int n = charToNum(c);
            n = r3.forwardTranslate(r2.forwardTranslate(r1.forwardTranslate(n))); //forward through
            n = reflector.forwardTranslate(n); //through reflector
            n = r1.backTranslate(r2.backTranslate(r3.backTranslate(n))); //backwards through
            r1.rotate();
            if((r1.turn%26) == r2.signal){
               r2.rotate();
            }
            if((r2.turn%26) == r3.signal){
               r3.rotate();
            }
            c = numToChar(n);
            output += c;         
         }
         
         System.out.println(output);
            
   }
   
   public static int charToNum(char c){
      return (int)c - 97;
   }
   
   public static char numToChar(int n){
      return (char)(n+97);
   }



}