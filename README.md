# Can s 3 y proteínas de transferencia de lípidos (nsLTP): co-sensibilización

**pMenor · Estadística clínica en R**

Código del análisis estadístico del artículo:

> Frías Ruiz C, Martín Martínez MD, García de la Rosa G. Cross-sensitization
> between *Cannabis sativa* and lipid transfer proteins: evidence from a large
> molecular cohort. *Allergy Asthma Proc.* 2026;47.
> [doi:10.2500/aap.2026.47.260045](https://doi.org/10.2500/aap.2026.47.260045)

## La pregunta

Can s 3 es la proteína de transferencia de lípidos (nsLTP) del cannabis. Por
su parecido estructural con las LTP de frutas, frutos secos, cereales y
pólenes, puede dar reacciones cruzadas, a veces graves (el llamado síndrome
cannabis-frutas-verduras). ¿Con qué frecuencia aparece la sensibilización a
Can s 3 y con qué otras nsLTP va asociada?

## El análisis

2.409 pacientes estudiados con la plataforma multiplex ALEX2: Can s 3, el
extracto de *C. sativa* y 12 nsLTP de alimentos y pólenes. Sensibilización:
IgE específica ≥ 0,35 kUA/L.

| Pregunta | Método |
|---|---|
| ¿Qué componentes predominan? | Prevalencia de sensibilización por componente |
| ¿Con qué va asociado Can s 3? | Tablas 2×2, test exacto de Fisher y odds ratio con IC 95 %, con ajuste de Bonferroni para 13 comparaciones |
| ¿Qué combinaciones son más frecuentes? | Gráfico UpSet |
| ¿Se mueven juntos los niveles de IgE? | Correlación de Spearman entre componentes |
| ¿Hay perfiles de pacientes? | Análisis de componentes principales y k-means (k = 3) |
| ¿Más componentes, más IgE total? | Correlación de Spearman |

Los valores fuera de los límites del laboratorio se convierten a número de
forma explícita: IgE específica "< 0,10" → 0,05 (la mitad del límite de
detección); IgE total "< 20" → 10 y "> 2500" → 2500. Solo afecta a los
análisis cuantitativos, no a quién se considera sensibilizado.

Para comprobar qué permitía la muestra se calculan la **precisión** de las
prevalencias y la **diferencia mínima detectable**, no la potencia a
posteriori, que depende del resultado obtenido.

### Resultados publicados

| Medida | Resultado |
|---|---|
| Sensibilización a Can s 3 | 3,7 % |
| nsLTP más frecuentes | Pru p 3 (11,7 %), Mal d 3 (10,0 %), Zea m 14 (8,8 %) |
| Asociaciones más fuertes con Can s 3 | Mal d 3 (OR 227,6) y Pru p 3 (OR 222,9) |
| Asociación más débil | Par j 2 (OR 2,41; no significativa tras Bonferroni) |
| Patrón | Can s 3 casi nunca aparece aislado: va dentro de una polisensibilización a nsLTP |

**Conclusión del artículo:** la sensibilización a Can s 3 forma parte de una
red de reactividad cruzada entre nsLTP, no es un fenómeno aislado. Por eso
conviene incluir Can s 3 en el estudio molecular de los pacientes con sospecha
de alergia a LTP.

## Datos

Este repositorio **no contiene datos de pacientes**. `R/simular_datos.R` genera
un dataset ficticio con la misma estructura (`data/cannabis_ltp_simulado.csv`),
con prevalencias y co-sensibilizaciones parecidas a las del estudio. Las
cifras se parecen a las publicadas, pero no son las mismas.

## Cómo reproducirlo

1. Clona el repositorio y abre la carpeta como proyecto de RStudio.
2. Instala los paquetes con las mismas versiones que se usaron (quedan
   registradas en `renv.lock`):

   ```r
   install.packages("renv")
   renv::restore()
   ```

3. Ejecuta `rmarkdown::render("analisis.Rmd")`.

El resultado, con código, salidas y figuras, se puede leer directamente en
[`analisis.md`](analisis.md). Las figuras usan las fuentes IBM Plex Sans e IBM
Plex Mono; si no están instaladas, R usa otras y avisa.

## Estructura

```
├── analisis.Rmd          análisis completo
├── analisis.md           el mismo análisis ya ejecutado (se lee en GitHub)
├── R/
│   ├── simular_datos.R   genera los datos ficticios
│   └── theme_pmenor.R    estilo de los gráficos
├── data/                 datos simulados
└── output/               figuras
```

## Licencia

El código se publica con licencia [MIT](LICENSE): puedes usarlo y adaptarlo citando la fuente. Si lo usas en un trabajo, cita también el artículo.

---

**pMenor** · Estadística clínica en R · [pmenor.es](https://pmenor.es)
