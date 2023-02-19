#include "io_text.h"

using namespace std;

string strip_line(string line) {
    string tmpstring1="";
    string tmpstring2="";
    unsigned short int flag=0,length;

    //write the line without preceding blanks
    for(unsigned int i=0; i<line.length(); i++) {
        if(line[i]!=' ') {
            flag=1;
        }
        if(flag==1) {
            tmpstring1+=line[i];
        }
    }

    //write the remaining line without concluding blanks  (inverted)
    flag=0;
    for(int i=tmpstring1.length()-1; i>=0; i--) {
        if(tmpstring1[i]!=' ') {
            flag=1;
        }
        if(flag==1) {
            tmpstring2+=tmpstring1[i];
        }
    }

    //reinvert the resulting string to get the correct order
    tmpstring1=tmpstring2;
    length=tmpstring2.length();
    for(int i=0; i<length; i++) {
        tmpstring1[length-i-1]=tmpstring2[i];
    }
    return tmpstring1;
}

string delete_char(string line,char del) {//self-explenatory
    string tmpstring="";
    for(unsigned int i=0; i<line.length(); i++) {
        if(line[i]!=del) {
            tmpstring+=line[i];
        }
    }
    return tmpstring;
}

string strip_blanks(string line) {
    //write only non-blank characters

    string tmpstring1="";
    for(unsigned int i=0; i<line.length(); i++) {
        if(line[i]!=' ') {
            tmpstring1+=line[i];
        }
    }
    return tmpstring1;
}
