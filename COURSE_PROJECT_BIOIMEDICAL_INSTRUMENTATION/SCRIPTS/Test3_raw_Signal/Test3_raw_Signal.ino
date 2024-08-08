float V_meas;
float V_internal = 4.99;
float R_pin = 38000.0;
float R_PSR;
long iTime_ave_usec = 1.0 / 60 * 1e6;
long iCounter;
long iTime0_usec;

void setup() {
  pinMode(A0, OUTPUT);
  pinMode(A1, OUTPUT);
  pinMode(A2, INPUT_PULLUP);
  pinMode(A3, OUTPUT);
  digitalWrite(A0, HIGH);
  digitalWrite(A1, LOW);
  digitalWrite(A3, LOW);
  Serial.begin(1000000);
}

void loop() {
  V_meas = 0.0;
  iCounter = 0;
  iTime0_usec = micros();
 
  while ((micros() - iTime0_usec) < iTime_ave_usec) {
    V_meas = analogRead(A2) / 1023.0 * 5.0 + V_meas;
    iCounter++;
  }
 
  V_meas = V_meas / iCounter;
  R_PSR = R_pin * V_meas / (V_internal - V_meas);
  Serial.print(R_PSR, 3);
  Serial.print(",");
  Serial.print(0);
  Serial.print(",");
  Serial.println(200000.0);
}
