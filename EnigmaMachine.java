import java.util.Scanner;
import javax.sound.sampled.Clip;

import java.io.File;
import java.io.ByteArrayInputStream;
import java.io.InputStream;
import java.io.IOException;

import java.net.URL;

import javax.sound.sampled.AudioFileFormat;
import javax.sound.sampled.AudioFormat;
import javax.sound.sampled.AudioInputStream;
import javax.sound.sampled.AudioSystem;
import javax.sound.sampled.DataLine;
import javax.sound.sampled.LineUnavailableException;
import javax.sound.sampled.SourceDataLine;
import javax.sound.sampled.UnsupportedAudioFileException;

import javax.sound.sampled.LineListener;
import javax.sound.sampled.LineEvent;

public class EnigmaMachine{

   public static void main(String args[]){
      /*Rotor r1 = new Rotor(1,3); //fast rotor
      Rotor r2 = new Rotor(2,2); //mid rotor
      Rotor r3 = new Rotor(3,5); //slow rotor
      */
      //reflector rotor
      StdAudio.loop("The_Imitation_Game.mid");
      
      Rotor reflector = new Rotor(0,1);
      reflector.rotor = new int[]{25,24,23,22,21,20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0};
      
      Scanner read = new Scanner(System.in);     
      
      //pick rotors and rotation settings
      System.out.println("Pick your 3 of 5 rotors (separated by a space):");
      String input = read.nextLine();    // get the entire line after the prompt 
      String[] rotors = input.split(" ");
      
      System.out.println("Pick your starting rotation for each rotors (separated by a space): ");
      input = read.nextLine();
      String[] rotation = input.split(" ");
      
      boolean escape = false;
            
      while(escape == false){
      
      Rotor r1 = new Rotor(Integer.valueOf(rotation[0]), Integer.valueOf(rotors[0]));
      Rotor r2 = new Rotor(Integer.valueOf(rotation[1]), Integer.valueOf(rotors[1]));
      Rotor r3 = new Rotor(Integer.valueOf(rotation[2]), Integer.valueOf(rotors[2]));
      
         //message prep
         System.out.println("\nPlease enter your message:");
         String test = read.nextLine();
         test = test.toLowerCase();
         test = test.replaceAll(" ","");
         test = test.replaceAll("\\p{P}","");
         //test = test.replace("\n", "");
         
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
            
   }
   
   public static int charToNum(char c){
      return (int)c - 97;
   }
   
   public static char numToChar(int n){
      return (char)(n+97);
   }

}