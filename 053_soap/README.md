# SOAP通信デモアプリケーション

PythonでSOAP通信を行うサーバーとクライアントのデモアプリケーションです。

## 必要なパッケージ

```bash
pip install -r requirements.txt
```

## 使い方

### 1. SOAPサーバーの起動

```bash
python soap_server.py
```

サーバーは `http://localhost:8000` で起動します。
WSDL は `http://localhost:8000/?wsdl` で確認できます。

### 2. SOAPクライアントの実行

別のターミナルで以下を実行:

```bash
python soap_client.py
```

## 提供されるサービス

### CalculatorService (計算サービス)

- `add(a, b)`: 2つの整数の足し算
- `subtract(a, b)`: 2つの整数の引き算
- `multiply(a, b)`: 2つの整数の掛け算
- `divide(a, b)`: 2つの数値の割り算

### UserService (ユーザー情報サービス)

- `get_user_name(user_id)`: ユーザー名を取得
- `get_user_email(user_id)`: メールアドレスを取得
- `get_user_age(user_id)`: 年齢を取得

## 使用ライブラリ

- **spyne**: SOAPサーバーの実装
- **zeep**: SOAPクライアントの実装
- **lxml**: XML処理
