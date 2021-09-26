/**
ガウスの消去法による連立一次方程式の数値計算法を行うDartライブラリ
A Dart Library for Solving Linear Equations by Gaussian Elimination
*/
/*
by Haruka Sato 佐藤陽花 (Takuya Matsunaga 松永拓也),
  日本大学理工学部応用情報工学科
  Division of Computer Engineering, College of Science and Technology, Nihon University.
  GitHub:@halka9000stg
  Mail:csta19097@g.nihon-u.ac.jp
*/
/**
@copyright 2021~ Haruka Sato under BSD3 License
*/


//整除関係
extension Divisibility on int{
  //整除
  int divisibility(int byNum){
    if(this.isDivisible(byNum)){
      return (this / byNum).floor();
    }else{
      return (this - (this % byNum)).divisibility(byNum);
    }
  }
  //整除可能か
  bool isDivisible(int byNum){
    return (this % byNum == 0);
  }
}
extension BasicOnInt on int{
  //符号反転
  int get reverse => -this;
}
extension FracCalc on int{
  Frac operator +(Frac y)=>Frac(this,1) + y;
  Frac operator -(Frac y)=>Frac(this, 1) - y;
  Frac operator *(Frac y)=>Frac(this, 1) * y;
  Frac operator /(Frac y)=>Frac(this, 1) / y;
}
class Frac{
  int _numer;
  int _denom;
  Frac(int numer, int denom){
    if(denom == 0){
      this._numer = 0;
      this._denom = -999;
      throw Exception("denominator is 0");
    }else if(numer == 0){
      this._numer = 0;
      this._denom = 1;
    }else{
      this._numer = numer;
      this._denom = denom;
    }
  }
int get numer => this._numer;
int get denom => this._denom;
  //約分
  Frac reduce(){
    int gcdNum = this._numer.gcd(this._denom);
    return this;
  }
  //通分
  Frac common(List<Frac> byFrac)=>[this,...byFrac].common()[0];
  //分母が0でないか
  Frac nonZeroCheck(){
    if(this._denom == 0){
      throw Exception("分母が0です");
    } else {
      return this;
    }
  }
  //整理：分母が0でないか確認する。その後分母が負の場合分母分子の符号をそれぞれ反転させる。さらに約分する。
  Frac regulate(){
    Frac nonZeroCheck = this.nonZeroCheck();
    if(nonZeroCheck._denom < 0){
      nonZeroCheck._numer = nonZeroCheck._numer.reverse;
      nonZeroCheck._denom = nonZeroCheck._denom.reverse;
    }
    return nonZeroCheck.reduce();
  }
  Frac operator +(Frac y){
    List<Frac> byFrac = [this, y];
    List<Frac> common = byFrac.common();
    return Frac(common[0]._numer + common[1]._numer, common[0]._denom).regulate();
  }
  Frac operator -(){
    this._numer.reverse;
    return this;
  }
  //逆数
  Frac operator ~(){
    return Frac(this._denom, this._numer).regulate();
  }
  Frac operator -(Frac y){
    return (this + (-y)).regulate();
  }
  Frac operator *(Frac y){
    return Frac(this._numer * y._numer, this._denom * y._denom).regulate();
  }
  Frac operator /(Frac y) => (this * (~y)).regulate();
  @override
  bool operator ==(Object other){
    if(other is Frac){
      return (this._numer == other._numer) && (this._denom == other._denom);
    }else{
      return false;
    }
  }
  bool operator <(Frac y){
    return (this._numer * y._denom) < (y._numer * this._denom);
  }
  bool operator >(Frac y){
    return (this._numer * y._denom) > (y._numer * this._denom);
  }
  bool operator <=(Frac y){
    return (this._numer * y._denom) <= (y._numer * this._denom);
  }
  bool operator >=(Frac y){
    return (this._numer * y._denom) >= (y._numer * this._denom);
  }
  @override
  String toString(){
    if(this._denom==1){
      return this._numer.toString();
    }else{
      return "${this._numer}/${this._denom}";
    }
  }
}
extension ListOfFrac on List<Frac>{
  //通分
  List<Frac> common(){
    List<int> commonDenom = this.map((Frac e)=>e.denom).toList();
    int commonDenomNum = commonDenom.reduce((a,b)=>a.lcm(b));
    List<Frac> commonFrac = this.map((Frac e)=>Frac((e.numer * commonDenomNum / e.denom).floor(), commonDenomNum)).toList();
    return commonFrac;
  }
}
class EquationLine{
  Frac _x;
  Frac _y;
  Frac _z;
  Frac _r;

  EquationLine(Frac x, Frac y, Frac z, Frac r){
    this._x = x;
    this._y = y;
    this._z = z;
    this._r = r;
  }

  Frac get x => this._x;
  Frac get y => this._y;
  Frac get z => this._z;
  Frac get r => this._r;

  EquationLine operator +(EquationLine y){
    return EquationLine(this.x + y.x, this.y + y.y, this.z + y.z, this.r + y.r);
  }
  EquationLine operator -(EquationLine y){
    return EquationLine(this.x - y.x, this.y - y.y, this.z - y.z, this.r - y.r);
  }
  EquationLine operator *(Frac y){
    return EquationLine(this.x * y, this.y * y, this.z * y, this.r * y);
  }
  EquationLine operator /(Frac y){
    return EquationLine(this.x / y, this.y / y, this.z / y, this.r / y);
  }
  @override
  String toString(){
    String hl = [this._x,this._y,this._z].map((Frac l)=>l.toString()).map((String l)=>"(\t$l)").join("+");
    String hr = "\t${this._r.toString()}";
    return [hl,hr].join("=");
  }
}
class Equation3{
  EquationLine _l1;
  EquationLine _l2;
  EquationLine _l3;

  Equation3(EquationLine l1, EquationLine l2, EquationLine l3){
    this._l1 = l1;
    this._l2 = l2;
    this._l3 = l3;
  }

  EquationLine get l1 => this._l1;
  EquationLine get l2 => this._l2;
  EquationLine get l3 => this._l3;

  Frac get pivot => this._l1.x;
  Equation3 pivotChoice(){
    if(this.pivot == Frac(0,1)){
      if(this._l2.x >= this._l3.x){
        this._changeL1L2();
      }else{
        this._changeL1L3();
      }
    }
    return this;
  }
  Equation3 _changeL1L2(){
    EquationLine temp = this._l1;
    this._l1 = this._l2;
    this._l2 = temp;
    return this;
  }
  Equation3 _changeL1L3(){
    EquationLine temp = this._l1;
    this._l1 = this._l3;
    this._l3 = temp;
    return this;
  }
  Map<String,Frac> solve(){
    this.pivotChoice();
    Map<String,Frac> result = {"x":Frac(0,1),"y":Frac(0,1),"z":Frac(0,1)};
    return result;
  }
  @override
  String toString(){
    return [this._l1,this._l2,this._l3].map((EquationLine l)=>l.toString()).join("\n");
  }
}