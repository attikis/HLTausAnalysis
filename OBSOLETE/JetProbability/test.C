void test()
{
  Int_t a = 3114;
  Int_t mod = a%10000;
  Int_t idd = 0;
  if(mod>1000)
    idd = mod/1000;
  else if (mod<1000)
    idd = mod/100;

  std::cout<<a<<"\t\t"<<mod<<"\t\t"<<idd<<endl;
}
