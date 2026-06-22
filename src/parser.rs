use logos::Logos;
use pyo3::prelude::*;
use pyo3_stub_gen::derive::gen_stub_pyfunction;

#[derive(Logos, Debug, PartialEq)]
#[logos(skip r"[ \t\n\f]+")]
enum Token {
    #[token("[")]
    LeftBracket,
    #[token("]")]
    RightBracket,
    #[regex(r#""[^"\n]*""#)]
    String,
    #[regex(r"<[^>\n]*>")]
    Path,
    #[regex(r#"[^ ;\]\[<>"\t\n\f]+"#)]
    Word,
    #[regex(r";[^\n]*", allow_greedy = true)]
    Comment,
}

#[gen_stub_pyfunction]
#[pyfunction]
/// Tokenize an input line.
///
/// * Will ignore comments starting with ;
/// * Tokens are separated by whitespace, or enclosed within <> or ""
/// * Additionally `[` and `]` are always separate tokens.
pub fn tokenize(line: String) -> Vec<(usize, usize)> {
    let mut res = vec![];
    let mut lex = Token::lexer(&line);
    while let Some(_) = lex.next() {
        res.push((lex.span().start, lex.span().end))
    }
    res
}
