[to English](https://github.com/rtomiyasu/ProjectPowderConograph/blob/main/README.md)

# Conograph CUI プログラムの操作説明書
以下では、オープンソースの粉末指数づけプログラム [Conogaph CUI Version 0.99](https://github.com/rtomiyasu/ProjectPowderConograph/tree/main/Conograph1_0_00_win)について簡単に説明します。Conograph GUI は現在開発中です (2013/3/2)。 Conograph は粉末指数づけの新しい数学的手法を採用しており、現時点でその手法を紹介している文献に、[[2](#References)]と[[3](#References)]があります。

<img alt="outline_JP" src="https://github.com/rtomiyasu/ProjectPowderConograph/assets/149344913/61c8335b-f1c1-4388-ae64-f13eff53a78e" width="40%">

```
図 1: Conograph による粉末指数づけの３つの主要ステージ
```
Conograph は中性子飛行時間法を含む任意の粉末回折データに対し、比較的短時間で解の徹底探索を行うことができます。 ピークサーチ実行後に行われる粉末指数づけは、図1の3つのステージに大別することができます。 Conograph は全てのブラベー格子、空間群、消滅則に共通する解の列挙手法を採用しているため、ブラベー格子決定が必要になります。 得られた格子定数は入力されたピーク位置の観測誤差から波及する誤差をある程度含むため、 誤差に安定なブラベー格子決定手法を新たに開発しました [[3](#References)]。

メモリ使用量優先探索 (regular search) を選択した場合、Conograph が準備している探索パラメータから特に大きな変更が行われなければ、全てのステージの実行は約10分程度で終了します。 速さ優先探索 (quick search) では、約5分以内に全ての処理が終了します。 (ただし、i7 CPU (3.2 GHz, 8スレッド)を使用した結果で、ご使用のパソコンによってはもっと時間がかかる場合があります。)

初めて Conograph を使用する際は、難しいケースでも入力する探索パラメータの変更の必要がほとんどない「メモリ使用量優先探索」を選択することをお勧めしています。ただし、「速さ優先探索」でも、ユニットセルが小さいまたは対称性が高いといった簡単なケースでは、 入力パラメータの変更が必要ないことが多いです。 「速さ優先探索」と「メモリ使用量優先探索」の違いは、[ここ](https://github.com/rtomiyasu/ProjectPowderConograph/blob/main/doc/TipsForPowderAutoIndexing_JP.md#Quick_searchとregular_searchの主な違い)でもう少し詳しく説明しています。

## NEWS
### 2016/9/7
- 底心格子に関する出力形式に関する誤りを訂正しました。

## FAQ
- [Conograph CUI の使い方](#Conograph_CUIの使い方)
    - [Conograph を実行する](#Conographを実行する)
    - [実行後に格子定数とゼロ点シフトを精密化する](#実行後に格子定数とゼロ点シフトを精密化する)
- [ピークサーチに関するアドバイス](#ピークサーチに関するアドバイス)
- [Conograph を用いた粉末指数づけについて](https://github.com/rtomiyasu/ProjectPowderConograph/blob/main/doc/TipsForPowderAutoIndexing_JP.md) (別ページに移動します)
    - 入力パラメータ表（XML ファイル）
    - 「速さ優先探索」と「メモリ使用量優先探索」の主な違い
    - 満足できる結果が得られないときに変更するパラメータ
    - 正しい解はどれですか? (figure of meritを用いた判定方法)
    - 粉末指数づけにおける解の一意性
- [CUI と GUI の主な違い](#CUIとGUIの主な違い)

### Conograph_CUIの使い方
#### Conographを実行する
1. 同プログラムを実行するには、以下の３つの入力ファイルを準備する必要があります。 (付属の "sample" フォルダに例があります。)
    - "*.inp.xml"：計算パラメータの入力ファイル ([例](https://github.com/rtomiyasu/ProjectPowderConograph/blob/main/Conograph1_0_00_win/sample/sample5/Cimetidine-SR.inp.xml)),
    - "cntl.inp.xml" ： "*.inp.xml"を含む入出力ファイル名を指定するファイル ([例](https://github.com/rtomiyasu/ProjectPowderConograph/blob/main/Conograph1_0_00_win/sample/sample5/cntl.inp.xml)),
    - IGOR テキストファイル: 粉末回折パターン (X, Y 座標と Y 座標の誤差) と以下のピークの情報を含むファイル (例。このファイルは付属の[ピークサーチプログラム](https://github.com/rtomiyasu/PeakSearch/tree/main)からも出力されます)
        1. ピーク位置 (2θ、time-of-flight、または d 値)、
        1. ピーク高さ (グラフ表示にのみ使用されます)、
        1. ピーク半値幅 (ピーク位置の誤差を推定するために使用します),
        1. 各ピークについて、指数づけおよび figure of merit の計算に使用しない/するを指定するための 0/1 フラグ。

1. "sample"フォルダの中の一つのフォルダをコピーし、 コピー先のフォルダ内の二つのXMLファイルの中身とファイル名"*.inp.xml"、IGOR テキストファイルの 0/1 フラグを適宜修正してください。 ファイル名"*.inp.xml"を変更した際は、"cntl.inp"のファイル名も併せて修正する必要があります。
1. コマンドプロンプト、またはお使いのOSのターミナルウィンドウを立ち上げ、 上で変更した"cntl.inp"と同じフォルダにカレントフォルダを移動してください。
1. "Conograph.exe"の絶対パスをコマンドラインから入力し、Conographを実行します。
1. CUI は格子定数のリストを含む XML ファイル（[六方晶の例](https://github.com/rtomiyasu/ProjectPowderConograph/blob/main/figures/sample2.index.xml)）を出力し、入力待ち状態に入ります。XML ファイルの最上部には、各々の figure of merit について最も良い値を得た格子定数がブラベー格子ごとに表示されます。そのすぐ下に、 得られた全ての格子定数の中で最も良い de Wolff figure of merit [[4](#References)] の値を得た格子定数が表示されます。

#### 実行後に格子定数とゼロ点シフトを精密化する
1. 出力 XML ファイルの中に表示されている格子定数は、"0403003"のような数字が付いています。 コマンドラインからその数字を入力すると、該当の格子定数とゼロ点シフトが精密化された後、以下の情報を含む IGOR テキストファイルが出力されます。([例](https://github.com/rtomiyasu/ProjectPowderConograph/blob/main/figures/sample2_lattice(Hexagonal_23.1%2C23.1%2C10.7%2C90%2C90%2C120_70.3).histogramIgor)):
    1. 入力 IGOR text file に含まれる情報のコピー、
    1. 指定した格子定数が与えるピーク位置、
    1. 指定した格子定数とブラベー格子。
1. 数字は何度でも入力できます。Conograph を終了させるにはquit" と入力してください。終了時に以下のファイルが出力されます:
    - XML ファイル: 精密化された各格子定数とゼロ点シフトの情報を記載したもの ([例](https://github.com/rtomiyasu/ProjectPowderConograph/blob/main/figures/sample2.index2.xml))
    - IGOR テキストファイル: 精密化の間に出力された IGOR テキストファイルを一つにまとめたもの ([例](https://github.com/rtomiyasu/ProjectPowderConograph/blob/main/figures/sample2_lattices.histogramIgor))

### ピークサーチに関するアドバイス
まずは「ピーク高さを基準に回折ピークを出来るだけ一様に拾うようにする」ことをお勧めしています。 （Conograph付属の[ピークサーチプログラム](https://github.com/rtomiyasu/PeakSearch/tree/main)を用いて、そのようなピークサーチ結果を得るためのパラメータの設定方法は付属の取扱説明書の中で紹介しています。） どの回折ピークを組み合わせればよいかはConographの列挙アルゴリズムによって比較的短時間の間に判定されます。根拠のある事前情報をお持ちであれば別ですが、そうではない場合、回折ピークの選別を行うことや重畳ピークを人為的に除外することで、入力情報の質を落とすことは避けた方がよいです。

### CUIとGUIの主な違い
ソフトウェア IGOR Pro をお持ちであれば、CUIを用いた解析もそれほど面倒ではありません。 以下では、GUIとの違いを説明します。

- GUI上ではピークサーチも実行できます。
- GUI 起動時に各パラメータの推奨値がテキストボックスに自動でセットされます。
- GUI起動時には recommended values are automatically set in the text boxes for respective input parameters.
- GUIでは、粉末回折パターンと各格子定数に対して計算されるピーク位置の比較を、より容易に高い自由度の下で行うことが出来ます。
- GUIでは、de Wolff figure of merit 以外の基準を用いて解のソートを行うことが出来ます。
- CUIでは一連の手続きとして実行している以下の機能を、 GUIでは独立に実行することが出来ますので、計算のやり直しに要するコストが小さくなります。>
    1. 指数づけ実行前に行う、reflection pair method [[1](#References)]によるゼロ点シフト Δ2θ の推定 (CUI使用時は、解析初期に画面上にいくつかの候補値が出力されます。)
    1. 解のソートに使用する de Wolff figure of merit Mnの計算に用いる n の推定 (dominant zone と呼ばれる問題が生じている場合、n = 20以上の値を使用する必要があります。)
- CUIでは、上記の推定値を使用する際、プログラムを一度終了させる必要があります。 特に a. について、 正しい値が Δ2θ = 0.195° と非常に大きなケースでも、 Conograph からは十分な結果が得られており、 多くのケースで Δ2θ = 0 を使用し、指数づけ後の精密化を行えば十分と考えられることから、CUIはこのような簡易な実装になっています。

図2 は CUI のフローチャートです:

<img alt="flowchart_JP" src="https://github.com/rtomiyasu/ProjectPowderConograph/assets/149344913/242c8e45-caf0-4d0b-b7db-2ddafe39cd60" width="40%">

```
図 2: CUIのフローチャート
```

## References
1. C. Dong, F. Wu, H. Chen,<br>Correction of zero shift in powder diffraction patterns using the reflection-pair method, J. Appl. Cryst., 32, pp. 850-853 (1999).
1. R. Oishi-Tomiyasu,<br>Distribution rules of crystallographic systematic absences on the Conway topograph and their application to powder auto-indexing, preprint.
1. R. Oishi-Tomiyasu,<br>Rapid Bravais-lattice determination algorithm for lattice constants containing large observation errors, Acta Cryst. A, 68, pp. 525-535 (2012).
1. P. M. de Wolff,<br>A simplified criterion for the reliability of a powder pattern indexing, J. Appl. Cryst., 1, pp. 108-113 (1968).
