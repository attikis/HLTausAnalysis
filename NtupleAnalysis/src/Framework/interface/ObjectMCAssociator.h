#ifndef ObjectMCAssociator_h
#define ObjectMCAssociator_h

class ObjectMCAssociator : public virtual MCTools, 
                           public virtual TreeDefinitionReco
{
 public: 
  short int FindTauMCProvenance(const unsigned short iTau);
  
  // private:
  
};

#endif //ObjectMCAssociator_h
