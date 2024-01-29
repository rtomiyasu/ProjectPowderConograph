[to English](https://github.com/rtomiyasu/ProjectPowderConograph/blob/main/doc/TipsForPowderAutoIndexing.md)

# Conographを用いた粉末指数づけについて
- [入力パラメータ表 (XML ファイル)](#入力パラメータ表_XML_ファイル)
- [Quick search と regular searchの主な違い](#Quick_searchとregular_searchの主な違い)
- [満足できる結果が得られないときに変更するパラメータ](#満足できる結果が得られないときに変更するパラメータ)
    1. より丁寧な探索を行う:
        - SearchLevel (XML),
        - CriticalValueForLinearSum (XML),
        - MaxNumberOfPeaks (XML).
    1. Figures of merit による解のソート結果を改善する:
        - ZeroPointShiftParameter (XML),
        - MaxNumberOfPeaksForFOM (XML),
        - 0/1 flags on use of the respective peaks for powder auto-indexing and computation of figures of merit (IGOR text).
    1. 計算を速くする:
        - NumberOfThreadsToUse (XML),
        - MinPrimitiveUnitCellVolume, MaxPrimitiveUnitCellVolume, MinDistanceBetweenLatticePoints (XML),
        - OutputCubicF, OutputCubicI, OutputCubicP, OutputHexagonal, OutputRhombohedral, OutputTetragonalI, OutputTetragonalP, - - OutputOrthorhombicF, OutputOrthorhombicI, OutputOrthorhombicC, OutputOrthorhombicP, OutputMonoclinicB, OutputMonoclinicP (XML),
        - SearchLevel & MaxNumberOfLatticeCandidates (XML).
- [正しい解はどれですか? (figure of meritを用いた判定方法)](#figure_of_meritを用いた判定方法)
- [粉末指数づけにおける解の一意性](#粉末指数づけにおける解の一意性)

## 入力パラメータ表_XML_ファイル
粉末指数づけを始める際は、まずは以下のテーブルにある推奨値（箱が灰色のものを除く）をそのまま使用してみてください。 パラメータ "AUTO" が入力されたとき、 推奨値は入力粉末回折パターンを用いて計算され、コマンドライン上に表示されます。
| XMLタグ | 説明 |推奨値 |
|----------|---------------|:-------------------:|
|回折計のパラメータ| | |
| IsAngleDispersion　　　　　| 0: 飛行時間法, 1: 角度分散法   |                   |
| ZeroPointShiftParameter　| 	(角度分散法のみ) ゼロ点シフト Δ2θ (°)  | 0 |
| WaveLength 　　　　　　　| (角度分散法のみ) 波長 (Å).　|                   |
| ConversionParameters    | (飛行時間法のみ) コンバージョンパラメータ、すなわち、定数項から始まる任意次数の多項式の系数 |                   |
|格子定数の列挙のためのパラメータ| | |
|SearchLevel|0: quick search (速さ優先探索。サイズが小さい、または、対称性の高いケースに向きます)、<br>1: regular search (メモリ使用量探索。難しい場合を含む全てのケースに対応します)。|    |
| MaxNumberOfPeaks | 使用するピークの最大数 | AUTO |
| CriticalValueForLinearSum | 伊藤の式を含むq (=1/d~^2)値の線形式が 0 に等しいかどうかを判定するための基準値 (i.e.｜ΣiQi｜<= cErr[ΣiQi]) | 1.0 |
| MinPrimitiveUnitCellVolume,<br>MaxPrimitiveUnitCellVolume | 	Primitive cell の体積の下限と上限 (Å3) | AUTO |
| MaxNumberOfTwoDimTopographs | 選別されたトポグラフから取得する伊藤の式を満たす4つのq値 (q1,q2,q3,q4) の最大数 | AUTO |
| MMaxNumberOfLatticeCandidates | 列挙する格子定数 (triclinic のみ) の最大数 | AUTO |
| ブラベー格子決定のためのパラメータ |  |  |
| OutputCubicF, OutputCubicI,<br>OutputCubicP, OutputHexagonal,<br>OutputRhombohedral, OutputTetragonalI,<br>OutputTetragonalP, OutputOrthorhombicF,<br>OutputOrthorhombicI, OutputOrthorhombicC,<br>OutputOrthorhombicP, OutputMonoclinicB,<br>OutputMonoclinicP, OutputTriclinic | 各々のブラベー格子について格子定数を出力する? (0:No, 1:Yes) | 1 |
| 格子定数のソートのためのパラメータ |  |  |
| Resolution | d* = 1/d値の相対的解像度<br>(二つの格子定数が与えるユニットセルがほぼ同一かどうかを判定するのに使用) | 0.03 |
| MaxNumberOfPeaksForFOM | figures of Merit の計算に使用するピークの数 n。 （回折パターンに含まれるピーク数が20より少なくても、20より小さい値に変更する必要はありません。） | 20 |
| MinFOM | de Wolff figure of merit Mn がこの値より大きい格子定数のみ出力されます。 | 3.0 |
| MinNumberOfMillerIndicesInRange,<br>MaxNumberOfMillerIndicesInRange | 各々の格子定数から、1本目から n 本目の観測されたピークの範囲に存在する、ピーク数を計算によって求めることが出来ます。そのピーク数の下限と上限。 | AUTO |
| MinUnitCellEdgeABC | 格子定数 a, b, c (Å)の下限 | 0 |
| MaxUnitCellEdgeABC | 格子定数 a, b, c (Å)の上限 | 1000 |
| 設定パラメータ |  | |
| NumberOfThreadsToUse | 	並列計算に使用するスレッド数 ("MAX" を入力すると、全スレッドが使用されます) | |
| AxisForRhombohedralSymmetry | "Rhombohedral" または "Hexagonal" を入力 (rhombohedral axis or hexagonal axis) | |
| AxisForMonoclinicSymmetry | "A" または "B" または "C" を入力 (A: a-axis, B : b-axis, C : c-axis). | B |
| ThresholdOnNormM | MnWu がこの値以下の格子定数はブラベー格子決定前に除去されます | 1.9 |
| ThresholdOnRevM | MnRev がこの値以下の格子定数はブラベー格子決定前に除去されます。 | 1.0 |
| MinDistanceBetweenLatticePoints | 結晶格子の最も距離が近い二つの格子点の距離 (Å)がこの値以下の格子定数は、格子定数の列挙を行うステージで除去されます。| 2.0 |

## Quick_searchとregular_searchの主な違い
二つの探索の違いは、 「計算速度とメモリ使用量のどちらを優先するか？」という技術的な問題から発生しており、 基礎となっている探索アルゴリズムは同じですので、 quick search において入力されたパラメータが適切であれば、二つの手法はほぼ同じ結果を返すと考えられます。 （ただし、quick search の方が2倍以上、計算時間が短くて済むでしょう。） 表2に二つの探索方法の違いをまとめます。

表 2: 二つの探索方法の違い

| | Quick search | Regular search |
|----------|------------|-------------|
| 計算時間 (Intel® Core™ i7 Processor (3.2 GHz), 8 スレッド使用) | < 5 分 | 約10分 |
| メモリ使用量 | 難しいケースでは、多くのメモリが使用されることがあり、 メモリ確保のエラーを予防するため、非常に多数の解が生成される場合、サイズの小さなユニットセルが優先的に列挙される格子定数の中から除去される | figure of merit の良い解から優先的に保存されるため、さほど多くは必要でない |
|入力パラメータの変更の必要性 |	小さな、または対称性の高いユニットセルではほとんど不要であるが、 それ以外のケースでは、以下のいずれかを変更する必要が生じる: MaxPrimitiveUnitCellVolume または MaxNumberOfLatticeCandidates | 多くのケースで不要 |

## 満足できる結果が得られないときに変更するパラメータ
1. より丁寧な探索を行う:
    - **SearchLevel: 0 → 1**<br>入力パラメータを推奨値にセットして Regular search を実行すれば、幅広いケースで十分に丁寧な探索を行うことができます。
    - **CriticalValueForLinearSum: 1 → 1.5**<br>特性X線または原子炉のデータにおいては、大き目の値を用いることで結果がよくなることがあるようです。
    - **MaxNumberOfPeaks: AUTO → 48よりも大きな数 (50--80)**<br>AUTOを指定した場合、48より大きな値が入ることはありません。 ピーク数を増やすと計算時間も増えますが、より多くの情報が粉末指数づけのために使われることになります。
1. Figures of merit による解のソート結果を改善する:
    - **ZeroPointShiftParameter: 0 → より正確なゼロ点シフト推定値**<br>ゼロ点シフトが大きな回折データでは (e.g., Δ2θ > 0.1 °) figures of merit の値がかなり小さくなることがあります。 粉末指数づけ後の格子定数とゼロ点シフトの精密化の結果が不十分な場合、 ゼロ点シフトの推定値を用いて粉末指数づけを行うことで、結果が改善されることがあります。 ゼロ点シフトの推定は以下の方法で行うことができます（ただしいずれの方法も間違えることがあるため、複数の候補を検討する必要があります）:
        1. (粉末指数づけ前) Reflection pair method を用いる [2]、
        1. (粉末指数づけ後) 正しそうな解または MnRev および MnSym が比較的大きな格子定数(例としてMnRev >3, MnSym >10程度) を用いて、ゼロ点シフトの精密化を行う。
    - **MaxNumberOfPeaksForFOM: 20 → 20よりも大きな数**<br>このパラメータには20がよく使われ、多くの標準的なケースで十分な結果を返します。 しかし、求めるユニットセルにdomnant zone と呼ばれるものが生じていると、20では不十分な場合があります。 格子定数とゼロ点シフトを精密化する際に、dominant zone が生じていれば、 コマンドラインで、警告メッセージと最低限必要なピーク数が表示されます。 ある格子定数を解だと判定する際は、コマンドライン上にこの警告メッセージが表示されていないことを確認してください。 もし表示されていた場合は、表示されているピーク数より大きな値をこのパラメータに対して入力する必要があります。
    - **ピークサーチからやり直す、または各ピークについて、指数づけおよび figure of merit の計算に使用しない/する、を指定するための 0/1 フラグ: 1 → 0**<br>不純物等が原因で生じる偽のピークは、 格子定数のソート結果に大きく影響することがあります。 （他方、格子定数の列挙は、不純物ピークの影響を受けにくいことが示されます。） 特に、パラメータ "MaxNumberOfPeaksForFOM" (figure of meritの計算に使用されるピーク数)が n のとき、 1本目から n 本目までのピークの間に不純物ピークが少なからず存在すると、 正しい解かどうかの判定は非常に困難になるため、 この範囲にあるピークについては、慎重にピークサーチを行う必要があります。 時には、ピークサーチからやり直すことや、 0/1フラグを用いて怪しいピークを除外することも有効です。 (ただし、ピークを除外しすぎると、格子定数の候補の列挙に影響が出てきます。)
1. 計算を速くする:
    - **NumberOfThreadsToUse**<br>一番単純な方法が、使用するスレッドの数を増やすことです。
    - **MinPrimitiveUnitCellVolume, MaxPrimitiveUnitCellVolume, MinDistanceBetweenLatticePoints**<br>上記のパラメータに関する情報を事前にお持ちであれば、格子定数の列挙にかかる時間を減らすことができます。 　(ブラベー格子ではなく、結晶の primitive cell (lattice) について、その体積、および最も距離が近い二つの格子点の距離を考える必要があります。)
    - **SearchLevel: 1 → 0 とし、さらに MaxNumberOfLatticeCandidates: AUTO → 64000 よりも大きな数 (100000 - 300000 程度)**<br>メモリ量の少ないコンピュータを考慮して、推奨値のAUTOによってセットされる MaxNumberOfLatticeCandidates の値は 64000 以下の数になります。 もし、この数よりも多数の格子定数が生成された場合、 体積の小さな格子定数から優先的に除去されるため、 quick search では、パラメータ MaxPrimitiveUnitCellVolume または MaxNumberOfLatticeCandidates の変更が必要になります。 しかし、例えば 4 GB のメモリを持つパソコンでは、MaxNumberOfLatticeCandidates = 300000 程度であれば、メモリ確保のエラーは起こらないようです。 このように、MaxNumberOfLatticeCandidates をご使用のパソコンに合わせて増やしておくことで、 quick search を用いて regular search と同等の結果を、 より多くのケースでより高速に得ることが出来ます。
    - **OutputCubicF, OutputCubicI, OutputCubicP, OutputHexagonal, OutputRhombohedral, OutputTetragonalI, OutputTetragonalP, OutputOrthorhombicF, OutputOrthorhombicI, OutputOrthorhombicC, OutputOrthorhombicP, OutputMonoclinicB, OutputMonoclinicP: 1 → 0**<br>もしブラベー格子について事前に情報をお持ちであれば、使用することで、ブラベー格子決定後のステージにかかる時間を減らすことができます。

## figure_of_meritを用いた判定方法
どの格子定数が正しいかを判定するための最も確実な方法は、画面上で観測ピーク位置と計算ピーク位置の比較を行うことですが、figure of merit を用いると、あらかじめ可能性が高い格子定数に候補を自動で絞っておくことが出来るため、正しい格子定数を見つけるまでに必要な時間を短縮することができます。

以下の5つが現在使用可能な figure of merit です:

- $`M_n`$ : de Wolff figure of merit [[1](#References)]、
- $`{M_n}^{W_u}`$ : Wu によって提案された改善版 figure of merit [[4](#References)]、
- $`{M_n}^{Rev}`$ : reversed figure of merit、
- $`{M_n}^{Sym}`$ : symmetric figure of merit、
- $`MN_ε`$ : 探索によって得られた解の中で、ほぼ同一なユニットセルの数。

ただし、n と ε は入力パラメータ MaxNumberOfPeaksForFOM、Resolution の値が入ります。 最初の4つの figure of merit は、観測ピーク位置と計算ピーク位置が相関を持たないとき1に近い値を取るように定義されています。 もし、$`M_n > 10`$, $`{M_n}^{W_u}>10`$, $`{M_n}^{Rev}>3`$ or $`{M_n}^{Sym}>30`$ のいずれかを満たす格子定数が見つかれば、正しい解を与えている可能性が高いでしょう。<br>上記の figures of merit のうち、de Wolff figure of merit $`M_n`$ が、常にではありませんが正しい解を検出する能力が最も高いです。 この理由として、de Wolff 評価指数の持つ以下の性質のうち、特に c が影響していると考えられます。

1. 観測されていない計算ピーク位置（消滅則などによる）の影響を受けにくい
1. 観測ピーク位置のうち指数がついてないもの（不純物ピークなどによる）の影響をうけやすい
1. ほぼ同一の格子定数が、異なるブラベー格子に含まれている場合、対称性の高いブラベー格子のものがより高い値を得る

$`{M_n}^{W_u}`$ も性質 a, b を持ちますが、c に関しては全く逆の傾向で、対称性の低いブラベー格子がより高い値を得ます。 $`{M_n}^{Rev}`$ の性質は、$`M_n`$ の上記の性質とは全て逆になり、特に不純物ピークの影響を受けにくいです。 $`{M_n}^{Sym}`$ は、$`M_n`$ と $`{M_n}^{Rev}`$ の中間の性質を持つように定義されています。

このことから、$`M_n`$が十分に大きく（例えば、$`M_n > 10`$ 程度）、 また他の複数の評価指数からも（少なくとも同じブラベー格子に属す解の中で）最上位に位置付けられている格子定数が存在すれば、 その格子定数は解である可能性はかなり高いと言えます。 同じブラベー格子に属す解の中で、異なる評価指数がそれぞれ別の格子定数を選んでいる場合、 多種多様な格子定数が生成されていると考えられるので、それぞれの評価指数がどのような格子定数を選んでいるのか、 一つ一つチェックした方がよいです。そのような解に対して画面上でピーク位置を確認することで最終判断ができる場合も多いでしょう。

ただし、仮に列挙された格子定数の中に正しい格子定数が存在し、 画面上でピーク位置の比較を行ったとしても、どれがその正しい解か区別できないことは起こり得ます。 よくある理由は粉末回折パターンの質が悪いということですが、 この問題については、[粉末指数づけにおける解の一意性](#粉末指数づけにおける解の一意性)もあわせて読まれることをお勧めします。

## 粉末指数づけにおける解の一意性
ピーク位置の情報からは、必ずしも格子定数が唯一つに決まるとは限りません。 複数の解が存在することは、低対称の格子定数では稀ですが、 立方晶、六方晶、三方晶の格子定数には、より対称性の低いブラベー格子に属する格子定数で、 同一のピーク位置を持つものが存在することが知られています [[3](#References)]。 これは数学によって示される現象のため、観測上の精度に関わりません。 以下の図は、cubic(F) の場合の具体例を示しています:

![ThreeSolutions](https://github.com/rtomiyasu/ProjectPowderConograph/assets/149344913/d5c9b94b-f412-4875-a841-c03a0f8fa74e)
```
図 1 : 異なる格子定数が完全に一致するピーク位置を持つ例
```
Conograph は徹底探索を行っているので、複数の解が存在する場合にも、 観測上の問題が大きい場合を除いて、全ての解が列挙された格子定数のリストに含まれていると考えられます。 その場合、最も対称性の高い解が、最良の de Wolff figure of merit の値を得ることが多いです。 （一般に対称性の高い解の方が、正しい確率が高いとされていますが、もちろん、常に正しいわけではありません。） 消滅則、または格子定数のサイズを考慮することで、格子定数の候補の数を減らせることがあります。

## References
1. P. M. de Wolff,<br>A simplified criterion for the reliability of a powder pattern indexing, J. Appl. Cryst., 1, pp. 108-113 (1968).
1. C. Dong, F. Wu, H. Chen,<br>Correction of zero shift in powder diffraction patterns using the reflection-pair method, J. Appl. Cryst., 32, pp. 850-853 (1999).
1. A. D. Mighell & A. Santoro,<br>Geometrical Ambiguities in the Indexing of Powder Patterns, J. Appl. Cryst., 8, pp. 372-374 (1975).
1. E. Wu,<br>A modification of the de Wolff figure of merit for reliability of powder pattern indexing, J. Appl. Cryst., 21, pp. 530-535 (1988).