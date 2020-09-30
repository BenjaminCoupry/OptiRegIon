using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using System.Reflection;
using System.Text.RegularExpressions;

namespace OptiRegIon
{
    class Program
    {
        static void Main(string[] args)
        {


            

            //Penser a faire un interpreteur de string (string -> Func<Vector,Vector,Vector>), puis a remplacer les params par les params optimaux trouves pour renvoyer la string
            Optimisateur op = new Optimisateur(0.00000001,0.1,1000,4,10,0.3,0.4);
            string fs = "p0 / (p1 + (Exp (p2 * x0))) ; Sin ( (p1 * x0) + p3 )";
            string fs2 = "p0 + (p3 * (Atan ((p1 * x0) + p2))) ; Sin ( (p4 * x0) + p5 )";
            string fs3 = "p0 + (p3 * ((Cos ((p1 * x0) + p2)) ^ 2 )) ; Sin ( (p4 * x0) + p5 )";
            /*
            List<Vector> x_ = new List<Vector>();
            for(int i=0;i<10;i++)
            {
                x_.Add(new Vector(new double[] { i }));
            }
            
            Console.WriteLine(Operation.CompterParams(fs));
            Func<Vector, Vector, Vector> F = Operation.Generer(fs);
            Vector paramreel = new Vector(new double[] { 0.7, 0.2, -0.6});
            Func<Vector, Vector> f = Optimisateur.Parametrer(F, paramreel);
            List<Tuple<Vector, Vector>> Attendu = Optimisateur.Evaluer(x_, f);
            */
            Console.WriteLine(op.FonctionOptimale(fs, "D:/lab/optimisation/att.csv"));
            Console.WriteLine(op.FonctionOptimale(fs2, "D:/lab/optimisation/att.csv"));
            Console.WriteLine(op.FonctionOptimale(fs3, "D:/lab/optimisation/att.csv"));
            /*
            Vector paramsaprox = op.ParametresOptimaux(F, Attendu, new Vector(new double[] { -1, -1,-1}), new Vector(new double[] { 1, 1,1}));
            Console.WriteLine(paramsaprox.ToString());
            Func<Vector, Vector> f_ = Optimisateur.Parametrer(F, paramsaprox);
            List<Tuple<Vector, Vector>> Obtenu = Optimisateur.Evaluer(x_, f_);
            Optimisateur.WriteToCsv(Attendu,"D:/lab/optimisation/att.csv");
            Optimisateur.WriteToCsv(Obtenu, "D:/lab/optimisation/obt.csv");
            */



        }
    }
    
    class Vector
    {
        int Dim;
        double[] valeurs;
        public Vector( int n)
        {
            Dim = n;
            valeurs = new double[n]; 
        }
        public Vector(double[] valeurs_)
        {
            valeurs = valeurs_;
            Dim = valeurs_.Length;
        }
        public static double Dist2(Vector A, Vector B)
        {
            return Math.Pow((A - B).Norme(), 2);
        }
        public static Vector operator + (Vector A, Vector B)
        {
            int n = A.GetDim();
            Vector ret = new Vector(n);
            for(int i=0;i<n;i++)
            {
                ret.Set(A.Get(i) + B.Get(i), i);
            }
            return ret;
        }
        public static double operator *(Vector A, Vector B)
        {
            int n = A.GetDim();
            double ret = 0;
            for (int i = 0; i < n; i++)
            {
                ret+= A.Get(i) * B.Get(i);
            }
            return ret;

        }
        public static Vector operator -(Vector A, Vector B)
        {
            int n = A.GetDim();
            Vector ret = new Vector(n);
            for (int i = 0; i < n; i++)
            {
                ret.Set(A.Get(i) - B.Get(i), i);
            }
            return ret;
        }
        public static Vector operator *(double k, Vector B)
        {
            int n = B.GetDim();
            Vector ret = new Vector(n);
            for (int i = 0; i < n; i++)
            {
                ret.Set(k * B.Get(i), i);
            }
            return ret;
        }
        public static Vector Unitaire(int N, int n)
        {
            Vector ret = new Vector(N);
            for (int i = 0; i < N; i++)
            {
                if (i == n)
                {
                    ret.Set(1,i);
                }
                else
                {
                    ret.Set(0, i);
                }
            }
            return ret;
        }
        public int GetDim()
        {
            return Dim;
        }
        public void Set(double[] Valeurs)
        {
            Dim = valeurs.Length;
            valeurs = Valeurs;
        }
        public void Set(double valeur, int i)
        {
            valeurs[i] = valeur;
        }
        public double Get(int i)
        {
            return valeurs[i];
        }
        public double Norme()
        {
            double val = 0;
            for(int i=0;i<Dim;i++)
            {
                val += Math.Pow(valeurs[i], 2);
            }
            return Math.Sqrt(val);
        }
        public Vector Normer()
        {
            double norme = this.Norme();
            if(norme ==0)
            {
                return Vector.Copie(this);
            }
            else
            {
                return (1.0 / norme) * this;
            }
            
        }
        public static Vector Rand(ref Random r, Vector min, Vector max)
        {
            int n = min.GetDim();
            Vector ret = new Vector(n);
            for (int i = 0; i < n; i++)
            {
                double val = min.Get(i) + r.NextDouble() * (max.Get(i) - min.Get(i));
                ret.Set(val, i);
            }
            return ret;
        }
        public static Vector Ones(int n)
        {
            
            Vector ret = new Vector(n);
            for (int i = 0; i < n; i++)
            {
                double val = 1;
                ret.Set(val, i);
            }
            return ret;
        }
        public static Vector Copie(Vector V)
        {
            int n = V.GetDim();
            Vector ret = new Vector(n);
            for (int i = 0; i < n; i++)
            {
                double val = V.Get(i);
                ret.Set(val, i);
            }
            return ret;
        }
        public static Vector FromString(string st)
        {
            st = st.Replace("(", "");
            st = st.Replace(")", "");
            string[] elt = st.Split(new char[] {' '});
            double[] res = elt.Select(x => Convert.ToDouble(x)).ToArray();
            Vector resv = new Vector(res.Count());
            resv.Set(res);
            return resv;
        }
        public override string ToString()
        {
            string result = "(";
            for(int i=0;i<Dim;i++)
            {
                double v = valeurs[i];
                
                result += v.ToString("F30");
                if(i!=Dim-1)
                {
                    result += " ";
                }
            }
            result += ")";
            return result;
        }
    }
    class Optimisateur
    {
        double delta;
        double Epsilon;
        int ItMaxGrad;
        int nbStagnations;
        int ItMaxAlea;
        double C1;
        double C2;
        Random r;

        public Optimisateur(double delta, double epsilon, int itMaxGrad, int nbStagnations, int itMaxAlea, double c1, double c2)
        {
            this.delta = delta;
            Epsilon = epsilon;
            ItMaxGrad = itMaxGrad;
            this.nbStagnations = nbStagnations;
            ItMaxAlea = itMaxAlea;
            C1 = c1;
            C2 = c2;
            r = new Random();
        }



        //Renvoie la fonction Gradient
        public Func<Vector, Vector> Gradient (Func<Vector, double> f)
        {
            return (x) =>
                {
                    int n = x.GetDim();
                    double f_ = f(x);
                    Vector ret = new Vector(n);
                    for(int i=0;i<n;i++)
                    {
                        double f_1 = f(x + delta * Vector.Unitaire(n, i));
                        double d_x = (f_1 - f_) / delta;
                        ret.Set(d_x, i);
                    }
                    return ret;
                };
        }
        //Conditions de Wolfe en partant de X dans la direction d
        public double TrouverPas(Vector X,Vector d, Func<Vector, double> f)
        {
            Func<Vector, Vector> Grad = Gradient(f);
            double a = 0;
            double t = 1;
            double b = double.PositiveInfinity;
            
            Vector GradLoc = Grad(X);
            double Floc = f(X);
            while (true)
            {
                
                double Fproj = f(X + t * d);
                Vector Gradproj = Grad(X + t * d);
                if(Math.Abs(b-a)<delta)
                {
                    break;
                }
                if(Fproj>Floc+C1*t*GradLoc*d)
                {
                    b = t;
                    t = 0.5 * (a + b);
                }
                else if(Gradproj*d<C2*GradLoc*d)
                {
                    a = t;
                    if(double.IsPositiveInfinity(b))
                    {
                        t = 2.0 * a;
                    }
                    else
                    {
                        t = 0.5 * (a + b);
                    }
                }
                else
                {
                    break;
                }
            }
            return t;
        }
        //Suit le gradient en paretant de X0
        public Vector DescenteGradient(Func<Vector, double> f,Vector X0)
        {
            Vector X = Vector.Copie(X0);
            Func<Vector, Vector> Grad = Gradient(f);
            int n = 0;
            Vector Grad_;
            double norgrad;
            do
            {
                Grad_ = Grad(X);
                norgrad = Grad_.Norme();
                Vector Dir = -1*Grad_.Normer();
                double pas = TrouverPas(X,Dir,f);
                X += pas * Dir;
                n++;
            } while (n < ItMaxGrad && norgrad>Epsilon) ;
            return X;
        }
        //Minimise la fonction f en prenant des x0 entre ParamTestMin et ParamTestMax
        public Vector Minimiser(Func<Vector,double> f, Vector ParamTestMin, Vector ParamTestMax)
        {
            Vector X0;
            Vector Xmin = new Vector(0);
            double valmin = double.PositiveInfinity;
            int nbEchec = 0;
            int n = 0;
            do
            {
                X0= Vector.Rand(ref r, ParamTestMin, ParamTestMax);
                Vector X = DescenteGradient(f, X0);
                double val = f(X);
                if(val<valmin)
                {
                    valmin = val;
                    Xmin = X;
                    nbEchec = 0;
                }
                else
                {
                    nbEchec++;
                }
                Console.WriteLine(n + "("+nbEchec+")");
                n++;
            }
            while (nbEchec < nbStagnations && n<ItMaxAlea);
            return Xmin;
        }
        //Renvoie la fonction f(x,param) avec ses parametres definis
        public static Func<Vector, Vector> Parametrer(Func<Vector,Vector,Vector> f, Vector param)
        {
            return (x) => f(x, param);
        }
        //Renvoie l'erreur de f(x) sur l'ensemble des couples (x,f(x)) attendus
        public static double CalculerErreur(List<Tuple<Vector,Vector>> Attendu, Func<Vector, Vector> f)
        {
            double sum = 0;
            foreach (Tuple<Vector, Vector> couple in Attendu)
            {
                sum += Vector.Dist2(f(couple.Item1), couple.Item2);
            }
            return sum;
        }
        //Donne l'erreur de f(x,param) sur l'ensemble des couples (x,f(x)) attendus en fonction des parametres
        public static Func<Vector,double> ErreurParametrique(Func<Vector, Vector, Vector> f, List<Tuple<Vector, Vector>> Attendu)
        {
            return (x) => CalculerErreur(Attendu, Parametrer(f, x));
        }
        // Renvoie les parametres qui minimisent l'erreur de f(x,param) sur l'ensemble des couples (x,f(x)) attendus
        public Vector ParametresOptimaux(Func<Vector, Vector, Vector> f, List<Tuple<Vector, Vector>> Attendu, Vector ParamTestMin, Vector ParamTestMax)
        {
            Func<Vector, double> errParam = ErreurParametrique(f, Attendu);
            Vector paramopti = Minimiser(errParam, ParamTestMin, ParamTestMax);
            Console.WriteLine("> Optimisation avec erreur de " + errParam(paramopti));
            return paramopti;
        }
        //Renvoie la fonction qui minimise l'erreur de f(x,param) sur l'ensemble des couples (x,f(x)) attendus
        public Func<Vector, Vector> FonctionOptimale(Func<Vector, Vector, Vector> f, List<Tuple<Vector, Vector>> Attendu, Vector ParamTestMin, Vector ParamTestMax)
        {
            return Parametrer(f, ParametresOptimaux(f, Attendu,ParamTestMin,ParamTestMax));
        }
        //Trouve l'expression de la fonction qui minimise l'erreur de f(x,param) sur l'ensemble des couples (x,f(x)) attendus
        public string FonctionOptimale(string f, List<Tuple<Vector, Vector>> Attendu, Vector ParamTestMin, Vector ParamTestMax)
        {
            string s = (string)f.Clone();
            Func<Vector, Vector, Vector>  f_eq = Operation.Generer(f);
            Vector paramsOP = ParametresOptimaux(f_eq, Attendu, ParamTestMin, ParamTestMax);
            for (int i = 0; i < paramsOP.GetDim(); i++)
            {
                s = s.Replace("p" + i, paramsOP.Get(i).ToString());
            }
            return s;

        }
        //Trouve l'expression de la fonction qui minimise l'erreur de f(x,param) sur l'ensemble des couples (x,f(x)) attendus
        public string FonctionOptimale(string f, List<Tuple<Vector, Vector>> Attendu)
        {
            int np = Operation.CompterParams(f);
            Vector ParamTestMin= -1*Vector.Ones(np);
            Vector ParamTestMax = Vector.Ones(np);
            string s = (string)f.Clone();
            Func<Vector, Vector, Vector> f_eq = Operation.Generer(f);
            Vector paramsOP = ParametresOptimaux(f_eq, Attendu, ParamTestMin, ParamTestMax);
            for (int i = 0; i < paramsOP.GetDim(); i++)
            {
                s = s.Replace("p" + i, paramsOP.Get(i).ToString());
            }
            return s;

        }
        //Trouve l'expression de la fonction qui minimise l'erreur de f(x,param) sur l'ensemble des couples (x,f(x)) stockes dans le csv specifie
        public string FonctionOptimale(string f, string PathAttendu)
        {
            return FonctionOptimale(f,ReadFromCsv(PathAttendu));
        }
        //Evalue une fonction en tous les points donnés par la liste abscisses
        public static List<Tuple<Vector, Vector>> Evaluer(List<Vector> abscisses, Func<Vector, Vector> f)
        {
            List<Tuple<Vector, Vector>> retour = new List<Tuple<Vector, Vector>>();
            foreach(Vector x in abscisses)
            {
                retour.Add(new Tuple<Vector, Vector>(x, f(x)));
            }
            return retour;
        }

        public static void WriteToCsv(List<Tuple<Vector, Vector>> Couples, string Path)
        {
            using (var fs = new StreamWriter(Path))
            {
                foreach(Tuple<Vector, Vector> T in Couples)
                {
                    fs.WriteLine(T.Item1.ToString()+";" +T.Item2.ToString());
                }
            }
        }
        public static List<Tuple<Vector, Vector>> ReadFromCsv(string Path)
        {
            List<Tuple<Vector, Vector>> result = new List<Tuple<Vector, Vector>>();
            string[] lines = System.IO.File.ReadAllLines(Path);
            foreach(string L in lines)
            {
                string[] u = L.Split(';');
                result.Add(new Tuple<Vector, Vector>(Vector.FromString(u[0]), Vector.FromString(u[1])));
            }
            return result;
        }
    }


    class Operation
    {
        static List<string> MathFunc = GetMathFunc();
        private static List<string> GetMathFunc()
        {
            List<string> ls = new List<string>();
            MethodInfo[] methodInfos = typeof(Math).GetMethods();
            foreach (MethodInfo mi in methodInfos)
            {
                if (!ls.Contains(mi.Name))
                {
                    ls.Add(mi.Name);
                }
            }
            return ls;
        }
        public static int CompterParams(string s)
        {
            MatchCollection count = new Regex("p[0-9]+").Matches(s);
            int max = -1;
            foreach(Match m in count)
            {
                string val = m.Value;
                int valp = Convert.ToInt32(val.Substring(1, val.Length - 1));
                if(valp>max)
                {
                    max = valp;
                }
            }
            return max+1;
        }
        public double EvaluerMathFunc(Dictionary<string, double> Variables)
        {
            if(MathFunc.Contains(Operateur))
            {
                Type[] typeOperandes = Operandes.Select(x => typeof(double)).ToArray();
                Object[] args = Operandes.Select(x => (Object)x.Evaluer(Variables)).ToArray();
                MethodInfo meth = typeof(Math).GetMethod(Operateur, typeOperandes);
                double resultat = (double)meth.Invoke(new Object(), args);
                return resultat;
            }
            else
            {
                return double.NaN;
            }
        }
        public bool isOp(string s)
        {
            return s == "+" || s=="-" || s == "*" || s == "^" || s == "/" || s=="=" || MathFunc.Contains(s);
        }
        public static List<string> TrouverBlocs(string s)
        {
            List<string> ret = new List<string>();
            bool DansParentheses = false;
            int nbParenth = 0;
            string curr = "";
            foreach (char c in s)
            {
                if (c == ')')
                {
                    nbParenth--;
                    if (nbParenth == 0)
                    {
                        DansParentheses = false;
                        if (curr != "" && curr != ")")
                        {
                            ret.Add(curr);
                        }
                        curr = "";
                    }
                }
                if (DansParentheses)
                {
                    curr += c;
                }
                else
                {
                    if(c != ' ' && c != '(')
                    {
                        curr += c;
                    }
                    if(c ==' ' || c=='(')
                    {
                        if (curr != "" && curr != ")")
                        {
                            ret.Add(curr);
                        }
                        curr = "";
                    }
                }
                if (c == '(')
                {
                    DansParentheses = true;
                    nbParenth++;
                }
                
            }
            if (curr != "" && curr != ")")
            {
                ret.Add(curr);
            }
            return ret;
        }
        public string Operateur;
        List<Operation> Operandes;
        public Operation (string s)
        {
            List<string> blocs = TrouverBlocs(s);
            Operandes = new List<Operation>();
            if(blocs.Count ==1 )
            {
                Operateur = blocs[0];
            }
            else
            {
                for (int i = 0; i < blocs.Count; i++)
                {
                    string sub = blocs.ElementAt(i);
                    if (sub != "")
                    {
                        if (isOp(sub))
                        {
                            Operateur = sub;
                        }
                        else
                        {
                            Operandes.Add(new Operation(sub));
                        }
                    }
                }
            }
            if(string.IsNullOrEmpty(s))
            {
                ;
            }

        }
        public string Ecrire(string tab = "")
        {
            string ret = "";
            ret += tab + "->"+Operateur + " : \n";
            for(int i=0;i<Operandes.Count;i++)
            {
                ret += Operandes[i].Ecrire(tab + "    ")+"\n";
            }
            return ret;
        }
        public double Evaluer()
        {
            return this.Evaluer(new Dictionary<string, double>());
        }
        public double Evaluer(Dictionary<string, double> Variables)
        {
            double valtmp;
            if(MathFunc.Contains(Operateur))
            {
                return EvaluerMathFunc(Variables);
            }
            if(double.TryParse(Operateur,out valtmp))
            {
                return valtmp;
            }
            else if(Variables.ContainsKey(Operateur))
            {
                return Variables[Operateur];
            }
            else
            {
                if(isOp(Operateur))
                {
                    double v;
                    switch(Operateur)
                    {
                        case "+":
                            v = Operandes.Sum((x) => x.Evaluer(Variables));
                            break;
                        case "-":
                            if(Operandes.Count == 1)
                            {
                                v = -Operandes[0].Evaluer(Variables);
                            }
                            else
                            {
                                v = Operandes[0].Evaluer(Variables) - Operandes[1].Evaluer(Variables);
                            }
                            break;
                        case "^":
                            v = Math.Pow(Operandes[0].Evaluer(Variables), Operandes[1].Evaluer(Variables));
                            break;
                        case "*":
                            v = Operandes[0].Evaluer(Variables) * Operandes[1].Evaluer(Variables);
                            break;
                        case "/":
                            v = Operandes[0].Evaluer(Variables) / Operandes[1].Evaluer(Variables);
                            break;
                        case "=":
                            v = Operandes[0].Evaluer(Variables);
                            break;

                        default:
                            v = double.NaN;
                            break;
                    }
                    return v;
                }
                else
                {
                    return double.NaN;
                }
            }
        }
        public void Remplacer(string old, string ns)
        {
            if(Operateur == old)
            {
                Operateur = ns;
            }
            for(int i=0;i<Operandes.Count;i++)
            {
                Operandes[i].Remplacer(old, ns);
            }
        }
        public static Func<Vector, Vector, Vector> Generer(string s)
        {
            Dictionary<string, double> Variables = new Dictionary<string, double>() ;
            string[] dimsres = s.Split(';');
            Operation[] ops = dimsres.Select(x_ => new Operation(x_)).ToArray();
            return (x, param) => {
                for (int i=0;i<x.GetDim();i++)
                {
                    Variables["x" + i] = x.Get(i);
                }
                for (int i = 0; i < param.GetDim(); i++)
                {
                    Variables["p" + i] = param.Get(i);
                }

                double[] values = ops.Select(x__ => x__.Evaluer(Variables)).ToArray();
                return new Vector(values);
            };
        }
    }
    
}
