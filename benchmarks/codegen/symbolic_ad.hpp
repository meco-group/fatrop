#pragma once
#include <cmath>
#include <iomanip>
#include <map>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>
#include <stdexcept>

// Symbolic graph for derivative generation.
namespace fatrop_benchmark::ad {
struct Symbol;
struct Graph {
    struct Node { std::string op; int a=-1, b=-1; double value=0; std::string name; };
    std::vector<Node> nodes;
    std::map<std::tuple<std::string,int,int,double,std::string>,int> unique;
    static inline Graph *active=nullptr;
    Graph() { active=this; }
    int intern(Node n) {
        auto key=std::make_tuple(n.op,n.a,n.b,n.value,n.name);
        auto found=unique.find(key);
        if(found!=unique.end()) return found->second;
        int i=int(nodes.size()); nodes.push_back(n); unique.emplace(key,i); return i;
    }
    int constant(double v) { return intern({"constant",-1,-1,v,{}}); }
    int input(const std::string &name) { return intern({"input",-1,-1,0,name}); }
    bool is(int a,double v) const { return nodes[a].op=="constant" && nodes[a].value==v; }
    int op(const std::string &op,int a,int b=-1) {
        if(op=="add") { if(is(a,0))return b; if(is(b,0))return a; }
        if(op=="sub") { if(is(b,0))return a; if(a==b)return constant(0); }
        if(op=="mul") { if(is(a,0)||is(b,0))return constant(0); if(is(a,1))return b; if(is(b,1))return a; }
        if(op=="div") { if(is(a,0))return constant(0); if(is(b,1))return a; }
        if(nodes[a].op=="constant" && (b<0 || nodes[b].op=="constant")) {
            double x=nodes[a].value, y=b<0?0:nodes[b].value;
            double v=op=="add"?x+y:op=="sub"?x-y:op=="mul"?x*y:op=="div"?x/y:
                op=="sin"?std::sin(x):op=="cos"?std::cos(x):op=="exp"?std::exp(x):op=="tan"?std::tan(x):std::tanh(x);
            if(std::isfinite(v))return constant(v);
        }
        return intern({op,a,b,0,{}});
    }
    std::pair<int,int> partials(int i) {
        const Node n=nodes[i]; int zero=constant(0),one=constant(1);
        if(n.op=="add")return {one,one};
        if(n.op=="sub")return {one,constant(-1)};
        if(n.op=="mul")return {n.b,n.a};
        if(n.op=="div")return {op("div",one,n.b),op("div",op("sub",zero,n.a),op("mul",n.b,n.b))};
        if(n.op=="sin")return {op("cos",n.a),zero};
        if(n.op=="cos")return {op("sub",zero,op("sin",n.a)),zero};
        if(n.op=="exp")return {i,zero};
        if(n.op=="tan")return {op("add",one,op("mul",i,i)),zero};
        if(n.op=="tanh")return {op("sub",one,op("mul",i,i)),zero};
        return {zero,zero};
    }
    std::vector<int> reverse(int output) {
        const int count=int(nodes.size()),zero=constant(0);
        std::vector<int> adj(count,zero);adj[output]=constant(1);
        for(int i=count-1;i>=0;--i) {
            Node n=nodes[i];if(n.a<0 || is(adj[i],0))continue;
            auto [da,db]=partials(i);
            adj[n.a]=op("add",adj[n.a],op("mul",adj[i],da));
            if(n.b>=0)adj[n.b]=op("add",adj[n.b],op("mul",adj[i],db));
        }
        return adj;
    }
    std::vector<int> forward(int variable,int count) {
        std::vector<int> dot(count,constant(0));dot[variable]=constant(1);
        for(int i=0;i<count;++i) {
            Node n=nodes[i];if(n.a<0)continue;
            auto [da,db]=partials(i);
            dot[i]=op("mul",da,dot[n.a]);
            if(n.b>=0)dot[i]=op("add",dot[i],op("mul",db,dot[n.b]));
        }
        return dot;
    }
    std::string emit(const std::vector<std::pair<std::string,int>> &outputs) const {
        std::vector<bool> used(nodes.size());
        auto mark=[&](auto &&self,int i)->void { if(used[i])return;used[i]=true;
            if(nodes[i].a>=0)self(self,nodes[i].a);if(nodes[i].b>=0)self(self,nodes[i].b); };
        for(auto &o:outputs)mark(mark,o.second);
        std::vector<std::string> refs(nodes.size());std::ostringstream code;int number=0;
        for(int i=0;i<int(nodes.size());++i)if(used[i]) {
            const auto &n=nodes[i];
            if(n.op=="constant") { std::ostringstream v;v<<std::setprecision(17)<<n.value;refs[i]=v.str();if(refs[i].find_first_of(".eE")==std::string::npos)refs[i]+=".0"; }
            else if(n.op=="input")refs[i]=n.name;
            else {
                refs[i]="v"+std::to_string(number++);std::string expr;
                if(n.b<0)expr="std::"+n.op+"("+refs[n.a]+")";
                else {std::string symbol=n.op=="add"?"+":n.op=="sub"?"-":n.op=="mul"?"*":"/";
                    expr="("+refs[n.a]+" "+symbol+" "+refs[n.b]+")";}
                code<<"    const double "<<refs[i]<<" = "<<expr<<";\n";
            }
        }
        for(auto &o:outputs)code<<"    "<<o.first<<" = "<<refs[o.second]<<";\n";
        return code.str();
    }
};
struct Symbol {
    int id;
    Symbol(double v=0):id(Graph::active->constant(v)){}
    static Symbol node(int i) { Symbol s;s.id=i;return s; }
};
inline Symbol operator+(Symbol a,Symbol b){return Symbol::node(Graph::active->op("add",a.id,b.id));}
inline Symbol operator-(Symbol a,Symbol b){return Symbol::node(Graph::active->op("sub",a.id,b.id));}
inline Symbol operator-(Symbol a){return Symbol(0)-a;}
inline Symbol operator*(Symbol a,Symbol b){return Symbol::node(Graph::active->op("mul",a.id,b.id));}
inline Symbol operator/(Symbol a,Symbol b){return Symbol::node(Graph::active->op("div",a.id,b.id));}
inline Symbol sin(Symbol a){return Symbol::node(Graph::active->op("sin",a.id));}
inline Symbol cos(Symbol a){return Symbol::node(Graph::active->op("cos",a.id));}
inline Symbol exp(Symbol a){return Symbol::node(Graph::active->op("exp",a.id));}
inline Symbol tan(Symbol a){return Symbol::node(Graph::active->op("tan",a.id));}
inline Symbol tanh(Symbol a){return Symbol::node(Graph::active->op("tanh",a.id));}
}
