#include <SoftwareSerial.h>
#include <TinyGPS.h>
#include <Wire.h>
int a1=0;
int a2=0;
int a3=0;
int a4=0;
float flat, flon;
unsigned long age;
int vout1;  
int x=0;
int in = 2;
int t1=0;
int cntt=0;
int count=0,i=0,k=0,rate=0;
unsigned long time2,time1;
unsigned long time;
SoftwareSerial ss(11,12);
TinyGPS gps;
void setup() 
{
ss.begin(9600);
Serial.begin(9600);
pinMode(10,INPUT); // Configuring pin A1 as input
pinMode(in, INPUT);
delay(1000);
}

void loop() 
{  
 t1=analogRead(A0);
 t1=t1/2.2;
 
vout1=analogRead(A1);

Serial.print(t1); 
Serial.print(",");   
Serial.print(vout1); 
Serial.print(",");   
Serial.println(digitalRead(10));

bool newData = false;
unsigned long chars;
unsigned short sentences, failed;
 k=0;
 cntt=cntt+1;
 Serial.println(cntt);
 while(k<5)
  {
  if(!digitalRead(in)) 
    {
      if(k==0)
      time1=millis();
      k++;
      while(!digitalRead(in));
    }
    }
    gps.f_get_position(&flat, &flon, &age);
    Serial.print("LAT=");
    Serial.print(flat == TinyGPS::GPS_INVALID_F_ANGLE ? 0.0 : flat, 6);
    Serial.print(" LON=");
    Serial.print(flon == TinyGPS::GPS_INVALID_F_ANGLE ? 0.0 : flon, 6);
 
      time2=millis();
      rate=time2-time1;
      rate=rate/5;
      rate=60000/rate;
      Serial.print(rate);
      Serial.print(",");
      Serial.println(t1);
  if(cntt>10)
  {
  cntt=0;
   ss.print("H:");
   ss.print(rate);
   ss.print(",");
   ss.print("T:");
   ss.println(t1);
   Serial.println("sent");
  }
if((vout1>400)&& (a1==0))
{   
ss.print("FALL DETECTED !");// turn the Buzzer off  
a1=1;
ss.print("T:");// turn the Buzzer off  
ss.print(t1);// turn the Buzzer off  
ss.print(", HB:");// turn the Buzzer off  
ss.print(rate);// turn the Buzzer off  
ss.print("  ");// turn the Buzzer off  
ss.print("http://maps.google.com/maps?q=");// turn the Buzzer off  
gps.f_get_position(&flat, &flon, &age);
ss.print(flat == TinyGPS::GPS_INVALID_F_ANGLE ? 0.0 : flat, 6);
ss.print(",");// The SMS text you want to send
ss.println(flon == TinyGPS::GPS_INVALID_F_ANGLE ? 0.0 : flon, 6);

}
if((t1>36) && (a2==0))
{   
ss.print("*high temperature ! ");// turn the Buzzer off  
a2=1;
ss.print("T:");// turn the Buzzer off  
ss.print(t1);// turn the Buzzer off  
ss.print(", HB:");// turn the Buzzer off  
ss.print(rate);// turn the Buzzer off  
ss.print("  ");// turn the Buzzer off  
ss.print("http://maps.google.com/maps?q=");// turn the Buzzer off  
gps.f_get_position(&flat, &flon, &age);
ss.print(flat == TinyGPS::GPS_INVALID_F_ANGLE ? 0.0 : flat, 6);
ss.print(",");// The SMS text you want to send
ss.println(flon == TinyGPS::GPS_INVALID_F_ANGLE ? 0.0 : flon, 6);

}
if((rate>100) && (a3==0))
{   
a3=1;

ss.print("high pulse detected !");// turn the Buzzer off  
ss.print("T:");// turn the Buzzer off  
ss.print(t1);// turn the Buzzer off  
ss.print(", HB:");// turn the Buzzer off  
ss.print(rate);// turn the Buzzer off  
ss.print("  ");// turn the Buzzer off  
ss.print("http://maps.google.com/maps?q=");// turn the Buzzer off  
gps.f_get_position(&flat, &flon, &age);
ss.print(flat == TinyGPS::GPS_INVALID_F_ANGLE ? 0.0 : flat, 6);
ss.print(",");// The SMS text you want to send
ss.println(flon == TinyGPS::GPS_INVALID_F_ANGLE ? 0.0 : flon, 6);
  
}
if((digitalRead(10)==1) && a4==0 )
{   
ss.print("need help!");// turn the Buzzer off  
a4=1;
ss.print("T:");// turn the Buzzer off  
ss.print(t1);// turn the Buzzer off  
ss.print(", HB:");// turn the Buzzer off  
ss.print(rate);// turn the Buzzer off  
ss.print("  ");// turn the Buzzer off  
ss.print("http://maps.google.com/maps?q=");// turn the Buzzer off  
gps.f_get_position(&flat, &flon, &age);
ss.print(flat == TinyGPS::GPS_INVALID_F_ANGLE ? 0.0 : flat, 6);
ss.print(",");// The SMS text you want to send
ss.println(flon == TinyGPS::GPS_INVALID_F_ANGLE ? 0.0 : flon, 6);

}
      k=0;
      rate=0;
for (unsigned long start = millis(); millis() - start < 1000;)
  {
    while (ss.available())
    {
      char c = ss.read();
      if (gps.encode(c)) // Did a new valid sentence come in?
        newData = true;
    }  
}

 if (newData)
  {
    gps.f_get_position(&flat, &flon, &age);
    Serial.print("LAT=");
    Serial.print(flat == TinyGPS::GPS_INVALID_F_ANGLE ? 0.0 : flat, 6);
    Serial.print(" LON=");
    Serial.print(flon == TinyGPS::GPS_INVALID_F_ANGLE ? 0.0 : flon, 6);
    }
   }
