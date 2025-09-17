#include <iostream>
#include <fstream>
void InputFiles() {
  std::ifstream infile("inputFileList.txt");
  string base = "outfiles";
  std::ofstream file1(base+"1/infile1.txt"),   file2(base+"2/infile2.txt"),   file3(base+"3/infile3.txt"),
                file4(base+"4/infile4.txt"),   file5(base+"5/infile5.txt"),   file6(base+"6/infile6.txt"),  
                file7(base+"7/infile7.txt"),   file8(base+"8/infile8.txt"),   file9(base+"9/infile9.txt"),
                file10(base+"10/infile10.txt"),file11(base+"11/infile11.txt"),file12(base+"12/infile12.txt"),
                file13(base+"13/infile13.txt"),file14(base+"14/infile14.txt"),file15(base+"15/infile15.txt"),
                file16(base+"16/infile16.txt"),file17(base+"17/infile17.txt"),file18(base+"18/infile18.txt"),
                file19(base+"19/infile19.txt"),file20(base+"20/infile20.txt"),file21(base+"21/infile21.txt");
  int count = 0;
  string line;
  while(getline(infile,line)) {
    count++;
    if(count < 101) file1<<line<<"\n";
    if(count > 100 && count < 201) file2<<line<<"\n";
    if(count > 200 && count < 301) file3<<line<<"\n";
    if(count > 300 && count < 401) file4<<line<<"\n";
    if(count > 400 && count < 501) file5<<line<<"\n";
    if(count > 500 && count < 601) file6<<line<<"\n";
    if(count > 600 && count < 701) file7<<line<<"\n";
    if(count > 700 && count < 801) file8<<line<<"\n";
    if(count > 800 && count < 901) file9<<line<<"\n";
    if(count > 900 && count < 1001) file10<<line<<"\n";
    if(count > 1000 && count < 1101) file11<<line<<"\n";
    if(count > 1100 && count < 1201) file12<<line<<"\n";
    if(count > 1200 && count < 1301) file13<<line<<"\n";
    if(count > 1300 && count < 1401) file14<<line<<"\n";
    if(count > 1400 && count < 1501) file15<<line<<"\n";
    if(count > 1500 && count < 1601) file16<<line<<"\n";
    if(count > 1600 && count < 1701) file17<<line<<"\n";
    if(count > 1700 && count < 1801) file18<<line<<"\n";
    if(count > 1800 && count < 1901) file19<<line<<"\n";
    if(count > 1900 && count < 2001) file20<<line<<"\n";
    if(count > 2000 && count < 2101) file21<<line<<"\n";
  }
}
