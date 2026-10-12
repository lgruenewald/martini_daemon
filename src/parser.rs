/*
Higher performance parser core.

Does the following:
1. reads files, split it into lines, merges lines with backslash
2. tokenizes each line, applying #define replacements, including ones that change the number of tokens
  * Will ignore comments starting with ;
  * Tokens are separated by whitespace, or enclosed within <> or ""
  * Additionally `[` and `]` are always separate tokens.
3. type checks the tokens in unwrap() calls
4. error message templating

martini_daemon/__parser/parser.py does the following
* handling of C preprocessor directives, itp sections/directives
* calling the right directive functions on tokenized lines
* resolving #include's to paths

The following Python-accessible classes are defined:
* Lexer
* Token
* TokenList
* TokenParseException

*/

use logos::Logos;
use pyo3::IntoPyObjectExt;
use pyo3::types::{PyAny, PyNone};
use pyo3::{
    exceptions::{PyException, PyValueError},
    prelude::*,
};
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pymethods};
use std::collections::HashMap;
use std::f64;

#[derive(Debug)]
#[gen_stub_pyclass]
#[pyclass(extends=PyException)]
pub struct TokenParseException {
    #[pyo3(get)]
    token: Option<Py<Token>>,
    #[pyo3(get)]
    message: String,
}

#[gen_stub_pymethods]
#[pymethods]
impl TokenParseException {
    #[new]
    pub fn new_err(token: Option<Bound<'_, Token>>, message: String) -> Self {
        let token = token.map(|x| x.unbind());
        Self { token, message }
    }
}

#[derive(Logos, Debug, PartialEq)]
#[logos(skip r"[ \t\n\f]+")]
#[pyclass]
pub enum TokenKind {
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

#[gen_stub_pyclass]
#[derive(Debug)]
#[pyclass]
pub struct Token {
    #[pyo3(get)]
    pub content: String,
    #[pyo3(get)]
    pub line_num: usize,
    #[pyo3(get)]
    pub path: String,
    #[pyo3(get)]
    pub start: usize,
    #[pyo3(get)]
    pub end: usize,
}

#[gen_stub_pyclass]
#[derive(Debug)]
#[pyclass(skip_from_py_object)]
pub struct TokenList {
    pub line: String,
    #[pyo3(get)]
    pub path: String,
    #[pyo3(get)]
    pub line_num: usize,
    // computed on first unwrap, since #define inbetween construction and unwrap
    // should still affect it
    pub tokens: Vec<Py<Token>>,
    #[pyo3(get)]
    pub comment: Option<String>,
}

#[gen_stub_pymethods]
#[pymethods]
impl TokenList {
    #[new]
    pub fn new(line: String, path: String, line_num: usize, py: Python<'_>) -> PyResult<TokenList> {
        let mut tokens = vec![];
        let mut lex = TokenKind::lexer(&line);
        let mut comment: Option<String> = None;
        while let Some(_) = lex.next() {
            let start = lex.span().start;
            let end = lex.span().end;
            let content = line[start..end].to_string();
            if !content.starts_with(";") {
                tokens.push(Py::new(
                    py,
                    Token {
                        content,
                        line_num,
                        path: path.clone(),
                        start,
                        end,
                    },
                )?);
            } else {
                comment = Some(content);
            }
        }

        Ok(TokenList {
            line,
            path,
            line_num,
            tokens,
            comment,
        })
    }

    pub fn __getitem__(&mut self, index: isize, py: Python<'_>) -> PyResult<Py<Token>> {
        let len = self.tokens.len() as isize;
        let index = match index {
            k if k < 0 && k > -len => (len + k) as usize,
            k if k < len => k as usize,
            _ => {
                return Err(PyValueError::new_err(format!(
                    "Token index to getitem is out of range. Index is {index}, len is {len}."
                )));
            }
        };
        Ok(self.tokens[index].clone_ref(py))
    }

    pub fn __len__(&mut self) -> usize {
        self.tokens.len()
    }

    pub fn get_line(&self) -> String {
        self.line.to_string()
    }

    #[pyo3(signature = (index, type_filter, mandatory=true, default=None))]
    pub fn unwrap<'py>(
        self_: Bound<'py, Self>,
        index: isize,
        type_filter: &str,
        mandatory: bool,
        default: Option<Bound<'py, PyAny>>,
    ) -> PyResult<Py<PyAny>> {
        if mandatory && default.is_some() {
            return Err(PyValueError::new_err(
                "Default value specified for mandatory unwrap.",
            ));
        }
        let py = self_.py();
        let self_ = self_.borrow();
        let len = self_.tokens.len() as isize;
        let index = match index {
            k if k < 0 && k > -len => (len - k) as usize,
            k if k >= 0 && k < len => k as usize,
            _ => {
                if !mandatory {
                    match default {
                        Some(x) => return x.into_py_any(py),
                        None => return PyNone::get(py).into_py_any(py),
                    };
                }

                let token = self_.tokens.last().map(|x| x.clone_ref(self_.py()));
                return Err(PyErr::new::<TokenParseException, _>((
                    token,
                    format!("Not enough tokens, expected token at index {}.", index),
                )));
            }
        };

        let token = &self_.tokens[index];
        let content = &token.borrow(self_.py()).content;

        return match type_filter {
            "int" => match content.parse::<i64>() {
                Ok(i) => i.into_py_any(py),
                _ => Err(PyErr::new::<TokenParseException, _>((
                    Some(token.clone_ref(py)),
                    format!("Expected number, got {}.", content),
                ))),
            },
            "index" => match content.parse::<u64>() {
                Ok(i) if i > 0 => (i - 1).into_py_any(py),
                _ => Err(PyErr::new::<TokenParseException, _>((
                    Some(token.clone_ref(py)),
                    format!("Expected index, got {}.", content),
                ))),
            },
            "float" => match content.parse::<f64>() {
                Ok(f) => f.into_py_any(py),
                _ => Err(PyErr::new::<TokenParseException, _>((
                    Some(token.clone_ref(py)),
                    format!("Expected float, got {}.", content),
                ))),
            },
            "positive" => match content.parse::<f64>() {
                Ok(f) if f > 0.0 => f.into_py_any(py),
                _ => Err(PyErr::new::<TokenParseException, _>((
                    Some(token.clone_ref(py)),
                    format!("Expected positive non-zero number, got {}.", content),
                ))),
            },
            "degree" => match content.parse::<f64>() {
                Ok(degrees) => {
                    let mut radians = f64::to_radians(degrees);
                    while radians < -f64::consts::PI {
                        radians += f64::consts::TAU;
                    }
                    while radians > f64::consts::PI {
                        radians -= f64::consts::TAU;
                    }
                    radians.into_py_any(py)
                }
                _ => Err(PyErr::new::<TokenParseException, _>((
                    Some(token.clone_ref(py)),
                    format!("Expected float (in degrees), got {}.", content),
                ))),
            },
            "word" | "raw" => content.into_py_any(py),
            "pattern" => content.replace("{", "[").replace("}", "]").into_py_any(py),
            "pair" => {
                let parts = content.split(":").collect::<Vec<_>>();
                if parts.len() != 2 {
                    return Err(PyErr::new::<TokenParseException, _>((
                        Some(token.clone_ref(py)),
                        format!("Expected a pair (index:name), got {}.", content),
                    )));
                }
                let Ok(index) = parts[0].parse::<usize>() else {
                    return Err(PyErr::new::<TokenParseException, _>((
                        Some(token.clone_ref(py)),
                        format!(
                            "Expected an index before ':', got {} (in {}).",
                            parts[0], content
                        ),
                    )));
                };
                (index - 1, parts[1]).into_py_any(py)
            }
            "string" => {
                if content.len() < 2
                    || !(content.starts_with('"') && content.ends_with('"')
                        || content.starts_with('<') && content.ends_with('>'))
                {
                    return Err(PyErr::new::<TokenParseException, _>((
                        Some(token.clone_ref(py)),
                        format!("Expected token of type string, got {}.", content),
                    )));
                }

                let len = content.len();
                content[1..len - 1].into_py_any(py)
            }
            _ => Err(PyValueError::new_err("Unknown type filter.")),
        };
    }

    pub fn assert_at_least(&self, i: usize, py: Python<'_>) -> PyResult<()> {
        if self.tokens.len() < i {
            return Err(PyErr::new::<TokenParseException, _>((
                self.tokens.last().map(|x| x.clone_ref(py)),
                format!(
                    "Line has too few ({}) tokens, expected at least {}.",
                    self.tokens.len(),
                    i
                ),
            )));
        }
        Ok(())
    }

    pub fn assert_no_more_than(&self, i: usize, py: Python<'_>) -> PyResult<()> {
        if self.tokens.len() > i {
            return Err(PyErr::new::<TokenParseException, _>((
                self.tokens.last().map(|x| x.clone_ref(py)),
                format!(
                    "Line has too many ({}) tokens, expected at most {}.",
                    self.tokens.len(),
                    i
                ),
            )));
        }
        Ok(())
    }
}

#[gen_stub_pyclass]
#[pyclass]
pub struct Lexer {
    pub defines: HashMap<String, String>,
}

#[gen_stub_pymethods]
#[pymethods]
impl Lexer {
    #[new]
    pub fn new() -> Self {
        Self {
            defines: HashMap::new(),
        }
    }

    pub fn tokenize_file(self_: Bound<'_, Lexer>, path: &str) -> PyResult<Vec<Py<TokenList>>> {
        // 1. reads files, split it into lines, merges lines with backslash
        let py = self_.py();

        let f = std::fs::read_to_string(path)?;
        let mut cumulative = "".to_string();
        let mut res = Vec::<Py<TokenList>>::new();
        for (linenum, line) in f.lines().enumerate() {
            let line = line.trim_start().trim_end();
            if line.len() == 0 {
                continue;
            }
            if let Some(no_backspace) = line.strip_suffix("\\") {
                cumulative.push_str(no_backspace);
                cumulative.push_str("\n");
            } else {
                cumulative.push_str(line);
                res.push(Py::new(
                    py,
                    TokenList::new(
                        std::mem::take(&mut cumulative),
                        path.to_string(),
                        linenum,
                        py,
                    )?,
                )?);

                cumulative = "".to_string();
            }
        }

        Ok(res)
    }

    pub fn add_define(&mut self, key: &str, value: &str) {
        self.defines.insert(key.to_string(), value.to_string());
    }

    pub fn get_define(&self, key: &str) -> Option<String> {
        self.defines.get(key).cloned()
    }

    pub fn rm_define(&mut self, key: &str) {
        self.defines.remove(key);
    }

    pub fn do_replacements(&self, token_list: Bound<'_, TokenList>) -> PyResult<()> {
        let py = token_list.py();
        let mut token_list = token_list.borrow_mut();

        let mut toks = std::mem::take(&mut token_list.tokens);

        // go token by token
        let mut i = 0;
        while i < toks.len() {
            let tok = toks[i].borrow(py);
            if let Some(val) = self.defines.get(&tok.content) {
                drop(tok);
                let old_tok = toks.remove(i);
                let old_tok_ref = old_tok.borrow(py);
                let line_num = old_tok_ref.line_num;
                let path = old_tok_ref.path.clone();
                let start = old_tok_ref.start;
                let end = old_tok_ref.end;
                let mut n = 0;
                for new_tok in val.split_ascii_whitespace() {
                    if new_tok.len() == 0 {
                        continue;
                    }
                    let element = Token {
                        content: new_tok.to_string(),
                        line_num,
                        path: path.clone(),
                        start,
                        end,
                    };
                    toks.insert(i + n, Py::new(py, element)?);
                    n += 1;
                }
            } else {
                i += 1;
            }
        }

        token_list.tokens = toks;

        Ok(())
    }
}
