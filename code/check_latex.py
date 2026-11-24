"""Static LaTeX source checks; no TeX engine or third-party package required."""
from collections import Counter
from pathlib import Path
import re
import sys

# TeX/LaTeX primitives and commands supplied by the declared packages.
CORE = set('chi kappa omega phi varrho mathfrak mathsf gamma theta tau sigma iota xi zeta dot ddot texorpdfstring cap cdots downarrow eta infty lim longrightarrow mapsto ne perp setminus sim uparrow '.split()) | set('documentclass usepackage newcommand renewcommand newtheorem begin end '
           'title author date maketitle thispagestyle pagestyle setlength '
           'parindent parskip emergencystretch textwidth thepage headrulewidth '
           'section subsection appendix label ref pageref cite bibitem '
           'noindent small normalsize centering caption emph textit textbf texttt '
           'LaTeX quad qquad left right bigl bigr Bigl Bigr mathbf frac binom sum prod int '
           'min max log det in subseteq ge le cdot ldots to longmapsto '
           'lfloor rfloor alpha beta delta ell lambda mu nu varepsilon bigotimes'.split())
PACKAGES = {
    'amsmath': set('numberwithin text operatorname eqref widetilde'.split()),
    'amssymb': set('mathbb'.split()),
    'latex': set('mathcal'.split()),
    'amsthm': set('theoremstyle'.split()),
    'xcolor': set('definecolor'.split()),
    'hyperref': set('hypersetup urlstyle url'.split()),
    'enumitem': set('setlist'.split()),
    'fancyhdr': set('fancyhf fancyhead fancyfoot'.split()),
    'booktabs': set('toprule midrule bottomrule'.split()),
    'longtable': set('endfirsthead endhead endfoot endlastfoot'.split()),
    'cleveref': set('cref Cref'.split()),
}
ENVIRONMENTS=set('document abstract theorem lemma proposition corollary conjecture '
                 'remark proof equation equation* align align* gather gather* '
                 'table tabular longtable center cases enumerate itemize thebibliography'.split())


def cleaned(source):
    # Preserve escaped percent signs and remove comments only at unescaped %.
    return re.sub(r'(?<!\\)%[^\n]*', '', source)


def inspect(source):
    errors=[]
    s=cleaned(source)
    brace=[]
    i=0
    while i<len(s):
        if s[i]=='\\':
            i+=1
            if i<len(s) and s[i].isalpha():
                while i<len(s) and s[i].isalpha(): i+=1
            else: i+=1
            continue
        if s[i]=='{': brace.append(i)
        if s[i]=='}':
            if not brace: errors.append('unmatched closing brace')
            else: brace.pop()
        i+=1
    if brace: errors.append('%d unmatched opening braces'%len(brace))
    preamble=s.split('\\begin{document}',1)[0]
    declared=set(re.findall(r'\\(?:newcommand|renewcommand)\s*\{\\([A-Za-z]+)\}',preamble))
    theorem_envs=set(re.findall(r'\\newtheorem\s*\{([^}]+)\}',preamble))
    active=set(['latex'])
    for group in re.findall(r'\\usepackage(?:\[[^]]*\])?\{([^}]+)\}',preamble):
        active.update(group.split(','))
    known=CORE|declared
    for package in active: known.update(PACKAGES.get(package,set()))
    unknown=sorted(set(re.findall(r'\\([A-Za-z]+)',s))-known)
    if unknown: errors.append('undefined control words: '+', '.join(unknown))
    stack=[]
    envs=[]
    for command,name in re.findall(r'\\(begin|end)\{([^}]+)\}',s):
        envs.append(name)
        if name not in ENVIRONMENTS|theorem_envs: errors.append('unknown environment '+name)
        if command=='begin': stack.append(name)
        elif not stack or stack[-1]!=name: errors.append('environment mismatch at '+name)
        else: stack.pop()
    if stack: errors.append('unclosed environments: '+', '.join(stack))
    labels=re.findall(r'\\label\{([^}]+)\}',s)
    bibitems=re.findall(r'\\bibitem(?:\[[^]]*\])?\{([^}]+)\}',s)
    refs=[]
    cites=[]
    for command,keys in re.findall(r'\\(ref|pageref|eqref|cref|Cref|cite)(?:\[[^]]*\])?\{([^}]+)\}',s):
        (cites if command=='cite' else refs).extend(k.strip() for k in keys.split(','))
    if set(refs)-set(labels): errors.append('missing labels: '+str(sorted(set(refs)-set(labels))))
    if set(cites)-set(bibitems): errors.append('missing bibliography keys: '+str(sorted(set(cites)-set(bibitems))))
    for name,items in [('labels',labels),('bibliography keys',bibitems)]:
        duplicate=[k for k,v in Counter(items).items() if v>1]
        if duplicate: errors.append('duplicate '+name+': '+str(duplicate))
    # Count inline/display TeX math delimiters outside escaped dollar signs.
    dollars=re.findall(r'(?<!\\)\$\$|(?<!\\)\$',s)
    if dollars.count('$')%2 or dollars.count('$$')%2: errors.append('unpaired dollar delimiter')
    for a,b in [('\\[','\\]'),('\\(','\\)')]:
        if len(re.findall(r'(?<!\\)'+re.escape(a),s)) != len(re.findall(r'(?<!\\)'+re.escape(b),s)): errors.append('unpaired display delimiter '+a)
    return errors,{'environments':len(envs)//2,'labels':len(labels),
                   'references':len(refs),'bibliography entries':len(bibitems),
                   'citations':len(cites),'control words':len(set(re.findall(r'\\([A-Za-z]+)',s)))}


def main():
    path=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parents[1]/'main.tex'
    s=path.read_text()
    errors,counts=inspect(s)
    # These mutations are in memory; they verify that the checker rejects each
    # required failure class without writing malformed source files.
    probes=[s+'{',s.replace('\\end{document}','\\end{table}',1),
            s.replace('\\end{document}','\\ref{missingLabel}\\end{document}',1),
            s.replace('\\end{document}','\\cite{missingEntry}\\end{document}',1),
            s.replace('\\end{document}','\\nonexistentMacro\\end{document}',1)]
    assert all(inspect(probe)[0] for probe in probes)
    if errors:
        print('\n'.join('FAIL: '+error for error in errors));sys.exit(1)
    print('PASS: balanced braces')
    print('PASS: matching begin/end environments (%d)'%counts['environments'])
    print('PASS: all references resolve (%d uses, %d labels)'%(counts['references'],counts['labels']))
    print('PASS: all citations resolve (%d uses, %d bibliography entries)'%(counts['citations'],counts['bibliography entries']))
    print('PASS: no undefined control words (%d distinct; standard/package commands and preamble definitions)'%counts['control words'])
    print('PASS: math delimiter counts and unique labels/bibliography keys')
    print('PASS: rejection probes for braces, environments, references, citations and macros (5)')
    print('LIMIT: static checking does not replace compilation or visual inspection.')

if __name__=='__main__':main()
