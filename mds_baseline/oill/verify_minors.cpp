// Independent exhaustive finite-field minor verification, without Sage.
#include <iostream>
#include <vector>
#include <stdexcept>
#include <cstdint>
using namespace std;
static int mul[256][256], inv[256];
static int multiply(int a,int b,int polynomial,int bits) {
    int r=0; while(b) {if(b&1)r^=a;b>>=1;a<<=1;if(a&(1<<bits))a^=polynomial;} return r;
}
static vector<vector<int>> combinations(int m,int k) {
    vector<vector<int>> result; vector<int> a(k);for(int i=0;i<k;i++)a[i]=i;
    while(true) {result.push_back(a);int p=k-1;while(p>=0 && a[p]==m-k+p)--p;if(p<0)break;++a[p];for(int j=p+1;j<k;j++)a[j]=a[j-1]+1;}
    return result;
}
int main() {
  try {
    int m,bits,poly;if(!(cin>>m>>bits>>poly)||m<1||m>12||bits<1||bits>8||poly<(1<<bits)||poly>=(1<<(bits+1))) throw runtime_error("Invalid parameters");
    int q=1<<bits;
    for(int a=0;a<q;a++)for(int b=0;b<q;b++)mul[a][b]=multiply(a,b,poly,bits);
    for(int a=1;a<q;a++) {for(int b=1;b<q;b++)if(mul[a][b]==1){inv[a]=b;break;}if(!inv[a])throw runtime_error("Polynomial does not define a field");}
    vector<vector<int>> matrix(m,vector<int>(m));for(auto&r:matrix)for(auto&v:r)if(!(cin>>v)||v<0||v>=q)throw runtime_error("Invalid matrix");
    uint64_t checked=0;
    for(int k=1;k<=m;k++) {
      auto comb=combinations(m,k);
      for(const auto&rows:comb)for(const auto&cols:comb) {
        int a[12][12];for(int i=0;i<k;i++)for(int j=0;j<k;j++)a[i][j]=matrix[rows[i]][cols[j]];
        bool singular=false;
        for(int c=0;c<k;c++) {
            int p=c;while(p<k && a[p][c]==0)++p;if(p==k){singular=true;break;}
            if(p!=c)for(int j=c;j<k;j++)swap(a[p][j],a[c][j]);
            for(int i=c+1;i<k;i++)if(a[i][c]) {
                int factor=mul[a[i][c]][inv[a[c][c]]];
                for(int j=c+1;j<k;j++)a[i][j]^=mul[factor][a[c][j]];
                a[i][c]=0;
            }
        }
        ++checked;
        if(singular) {cout<<"{\"mds\":false,\"checked\":"<<checked<<",\"minor_size\":"<<k<<"}\n";return 0;}
      }
    }
    cout<<"{\"mds\":true,\"checked\":"<<checked<<"}\n";
    return 0;
  }catch(const exception&e){cerr<<e.what()<<"\n";return 2;}
}
