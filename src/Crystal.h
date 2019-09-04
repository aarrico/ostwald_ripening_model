class Crystal{
  public:
    Crystal(const char *);
    Crystal& operator=(const Crystal&);
    ~Crystal();
    void copyCrystal(Crystal, int);
    void fcn(double, double *);
    int test(double *);
    void createParams(const char *, double);
  private:  
    int nCrystals;  
    double gamma, c0, cStar, k, *mu, *xStar;
};

