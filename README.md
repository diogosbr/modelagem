# Modelagem de nicho
Funções para auxiliar a modelagem de distribuição de nicho

## **Funções** <p>
### *cut.raster*
Função para recortar as variáveis ambientais

### *cor.data*
Função para verificar a correlação entre as variáveis ambientais

### *clean*
Função para selecionar os pontos espacialmente únicos e retirar os NAs

### *modelos*
Função que roda os algoritmos

---
# Como usar
###cut.raster (raster.dir,shape.dir,extension=".asc",plot=F,trim=F)

Argumentos:

 * raster.dir: diretório que contém os rasters a serem cortados. Se não for informado, então vai procurar os rasters na pasta local.
 * shape.dir: diretório que stá o shape que será usado como máscara para cortar os rasters. Obrigatório.
 * extension: conjuntos de caracteres com a extensão dos ratsers. ".asc" é o padrão.
 * plot: plota os rasters. TRUE é o padrão.
 * trim: se for TRUE os NAs gerados após o corte dos rasters são removidos. FALSE é o padrão.

Exemplo de uso:

I'm no good at writing sample / filler text, so go write something yourself.

Look, a list!

 * foo
 * bar
 * baz

And here's some code! :+1:

```javascript
$(function(){
  $('div').html('I am a div.');
});
```

This is [on GitHub](https://github.com/jbt/markdown-editor) so let me know if I've b0rked it somewhere.


Props to Mr. Doob and his [code editor](http://mrdoob.com/projects/code-editor/), from which
the inspiration to this, and some handy implementation hints, came.

### Stuff used to make this:

 * [markdown-it](https://github.com/markdown-it/markdown-it) for Markdown parsing
 * [CodeMirror](http://codemirror.net/) for the awesome syntax-highlighted editor
 * [highlight.js](http://softwaremaniacs.org/soft/highlight/en/) for syntax highlighting in output code blocks
 * [js-deflate](https://github.com/dankogai/js-deflate) for gzipping of data to make it fit in URLs
