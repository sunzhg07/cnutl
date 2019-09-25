#ifndef _PROGRESSBAR_
#define _PROGRESSBAR_

#include <iostream>

// Packed by zhuo. 2017-11-08
// Usage 1> init(point_number, print_scale)
//          initial this class before your loop
//       2> PrintBar( point_pos ) point_pos = loop index.
//          Use it to print progress bar in your loop
//       3> EndBar cout<<endl; after your loop.

class Pgs_bar{
 public:
  int pos,pos_p,plan_t,num,scale;
  void init(int point_num, int print_scale){
    plan_t=0;
    pos = 0;
    scale=print_scale;
    num = point_num;
  };
  void PrintBar(int point_pos){
    if(point_pos > num)
      std::cout<<"point_pos Larger than point_num @ Pgs_bar::PringBar"<<endl;
    plan_t ++;
    pos = int( (point_pos+1)*scale/(1.0*num) );
    if( (plan_t*100)/(1.0*num) >2 ){
      std::cout<<"process [";
      for(int k=0; k<scale;k++){
	if(k<pos)cout<<"==";
	else if (k == pos) cout<<">";
	else std::cout << "  ";	
      }
      //plan_t = 0;
      std::cout << "] " << int( (point_pos+1)*100.0 /num)  << " % \r";
      std::cout.flush();
    }
  };
  void EndBar(){
    std::cout<<endl;
  };
};

#endif
