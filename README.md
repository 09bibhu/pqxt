# pqxt
Prototype Implementation for PQXT

//Compile and Run
1. Run "Make"
2. For PQXT, run ./pqxt
3. For OQXT, run ./oqxt_eurosnp

//Parameters in the uploaded version:
1. The number of keywords is 5
2. The maximum frequency for any keyword is 10000.

//Updating number of keywords and maximum frequency
1. Navigate to tools.h
2. Update "sizew" for changing the number of keywords
3. Update "sized" for changing the maximum frequency of any keyword
4. Update sizewid to (sizew times sizeid)

//Note
1. You need OpenSSL version 3.0.2 or more installed in your system
